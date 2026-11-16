"""Corridor roads as driven: each road routed through the OpenStreetMap network between its towns, in order.

    python src/01_roads.py

Reads data/raw/osm_corridor.gpkg (ogr2ogr from the Geofabrik Pakistan extract, box in config/corridor.toml).
For each road in config/corridor.toml: build a graph from the OSM ways of the allowed classes (nodes are shared OSM
vertices, so junctions connect; edges weighted by length), snap each town to the nearest graph node, and join the
shortest paths town to town. This replaces selecting by ref or name, which picked unrelated roads that share a name.
Each road also keeps which OSM ways it used, so bridges and tunnels come with it, and every 100 m segment records the
stretch it is on ("Abbottabad - Mansehra"). Snap distances are reported: a town more than 1 km from the network is
flagged, because its coordinate or the allowed classes are wrong.
Writes data/processed/corridor.gpkg (layers roads, segments, structures, dams, towns), outputs/roads_summary.csv.
"""
import re
import sys
import tomllib
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import linemerge, substring

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
CFG = tomllib.loads((R / "config/corridor.toml").read_text(encoding="utf-8"))
UTM = 32643


def tag(s, k):
    m = re.search(rf'"{k}"=>"([^"]*)"', s or "")
    return m.group(1) if m else None


def load_places():
    p = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="places").to_crs(UTM)
    p["en"] = p.other_tags.fillna("").str.extract(r'"name:en"=>"([^"]*)"')[0].fillna(p["name"]).fillna("")
    p["en"] = p["en"].str.replace("‎", "", regex=False).str.strip().str.lower()
    p["nm"] = p["name"].fillna("").str.strip().str.lower()
    return p


def geocode(v, places, prev):
    """A place from OSM by exact English or local name; nearest to the previous place if repeated; hint [lon, lat]
    preferred when given and an OSM match lies within 8 km of it (else the hint itself)."""
    name = v[0].split(" (")[0].strip().lower()
    m = places[(places.en == name) | (places.nm == name)]
    if len(v) == 3:
        h = gpd.GeoSeries(gpd.points_from_xy([v[1]], [v[2]]), crs=4326).to_crs(UTM).iloc[0]
        if len(m):
            d = m.distance(h)
            if d.min() < 8000:
                return m.geometry.loc[d.idxmin()]
        return h
    if not len(m):
        return None
    if prev is not None and len(m) > 1:
        return m.geometry.loc[m.distance(prev).idxmin()]
    return m.geometry.iloc[0]


def graph(ways):
    G = nx.Graph()
    for idx, g in zip(ways.index, ways.geometry):
        lines = [g] if g.geom_type == "LineString" else list(g.geoms)
        for ln in lines:
            c = [(round(x, 2), round(y, 2)) for x, y in ln.coords]
            for a, b in zip(c, c[1:]):
                if a != b:
                    G.add_edge(a, b, w=float(np.hypot(a[0] - b[0], a[1] - b[1])), way=idx)
    # OSM in remote valleys is often drawn in pieces whose ends do not quite meet: join a dead end to the nearest node
    # of another piece within 40 m (logged as a gap edge, weight x 3 so real links are preferred)
    from scipy.spatial import cKDTree
    nodes = list(G.nodes); xy = np.array(nodes); tree = cKDTree(xy)
    comp = {n: i for i, c in enumerate(nx.connected_components(G)) for n in c}
    gaps = 0
    for n in [n for n in nodes if G.degree(n) == 1]:
        for j in tree.query_ball_point(n, 40):
            m = nodes[j]
            if comp[m] != comp[n]:
                G.add_edge(n, m, w=3 * float(np.hypot(n[0] - m[0], n[1] - m[1])), way=None); gaps += 1
                break
    G.graph["gaps"] = gaps
    return G


def main():
    r = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="roads").to_crs(UTM)
    ot = r.other_tags.fillna("")
    r["ref"] = ot.map(lambda s: tag(s, "ref"))
    r["bridge"] = ot.str.contains('"bridge"=>"yes"|"bridge"=>"viaduct"')
    r["tunnel"] = ot.str.contains('"tunnel"=>"yes"')
    r["name_en"] = ot.map(lambda s: tag(s, "name:en")).fillna(r["name"])
    places = load_places()
    GA = None
    roads, segs, towns, summ = [], [], [], []
    for rd in CFG["road"]:
        refs = {x.replace("-", "").upper() for x in rd.get("refs", [])}
        by_ref = r.ref.fillna("").str.replace("-", "").str.upper().str.split(";").map(lambda v: bool(refs & set(v)))
        base = set(rd["classes"]) | ({"trunk", "trunk_link", "primary", "primary_link", "secondary"} if rd.get("group") != "trunk" else set())
        ways = r[r.highway.isin(base) | by_ref]          # allowed classes (+ main roads, so valleys connect), + the road's number
        G = graph(ways)
        nodes = np.array(list(G.nodes))
        pts = [geocode(v, places, prev) for prev, v in zip([None] + [None] * (len(rd["via"]) - 1), rd["via"])]
        for i in range(1, len(pts)):                     # re-resolve name-only places nearest to the previous one
            if len(rd["via"][i]) == 1:
                pts[i] = geocode(rd["via"][i], places, pts[i - 1])
        snapped = []
        for v, p in zip(rd["via"], pts):
            nm = v[0]
            if p is None:
                print(f"  ! {rd['id']}: {nm} not found in OSM places and no hint given", flush=True); continue
            d = np.hypot(nodes[:, 0] - p.x, nodes[:, 1] - p.y); i = int(d.argmin())
            snapped.append(tuple(nodes[i]))
            towns.append({"road": rd["id"], "town": nm, "snap_m": round(float(d[i])), "geometry": Point(nodes[i])})
            if d[i] > 1000:
                print(f"  ! {rd['id']}: {nm} is {d[i] / 1000:.1f} km from the {'/'.join(rd['classes'][:2])} network", flush=True)
        used = []
        names = [v[0] for v, p in zip(rd["via"], pts) if p is not None]
        for k, (a, b) in enumerate(zip(snapped, snapped[1:])):
            if a == b:
                continue
            try:
                path = nx.shortest_path(G, a, b, weight="w")
            except nx.NetworkXNoPath:
                if GA is None:
                    GA = graph(r)                                  # every road class in the extract
                    print(f"  all-class graph: {GA.number_of_nodes():,} nodes, {GA.graph['gaps']} gap joins", flush=True)
                na = np.array(list(GA.nodes))
                aa = tuple(na[np.hypot(na[:, 0] - a[0], na[:, 1] - a[1]).argmin()]); bb = tuple(na[np.hypot(na[:, 0] - b[0], na[:, 1] - b[1]).argmin()])
                try:
                    path = nx.shortest_path(GA, aa, bb, weight="w"); Gp = GA
                    print(f"  {rd['id']}: {names[k]} -> {names[k + 1]} routed over all road classes", flush=True)
                except nx.NetworkXNoPath:
                    print(f"  ! {rd['id']}: no path {names[k]} -> {names[k + 1]} even over all classes", flush=True); continue
            else:
                Gp = G
            if False:
                print(f"  ! {rd['id']}: no path {names[k]} -> {names[k + 1]}", flush=True); continue
            line = LineString(path)
            stretch = f"{names[k]} - {names[k + 1]}"
            ways_used = {Gp[u][v]["way"] for u, v in zip(path, path[1:])} - {None}
            used += list(ways_used)
            br = r.loc[list(ways_used)]
            roads.append({"road": rd["id"], "group": rd.get("group", ""), "road_name": rd["name"], "stretch": stretch, "km": line.length / 1000,
                          "bridge_km": br[br.bridge].length.sum() / 1000, "tunnel_km": br[br.tunnel].length.sum() / 1000, "geometry": line})
            n = max(1, int(np.ceil(line.length / 100)))
            sub = br[br.bridge | br.tunnel]
            for i in range(n):
                s = substring(line, i * line.length / n, (i + 1) * line.length / n)
                hit = sub[sub.distance(s.interpolate(0.5, normalized=True)) < 15] if len(sub) else sub
                segs.append({"road": rd["id"], "group": rd.get("group", ""), "stretch": stretch, "bridge": bool(len(hit) and hit.bridge.any()),
                             "tunnel": bool(len(hit) and hit.tunnel.any()), "geometry": s})
        rr = [x for x in roads if x["road"] == rd["id"]]
        summ.append({"road": rd["name"], "group": rd.get("group", ""), "km": round(sum(x["km"] for x in rr), 1), "stretches": len(rr),
                     "bridge_km": round(sum(x["bridge_km"] for x in rr), 2), "tunnel_km": round(sum(x["tunnel_km"] for x in rr), 2),
                     "worst_snap_m": max([t["snap_m"] for t in towns if t["road"] == rd["id"]] or [None]),
                     "worst_snap_town": max((t for t in towns if t["road"] == rd["id"]), key=lambda t: t["snap_m"], default={"town": None})["town"]})
        print(summ[-1], flush=True)
    roads = gpd.GeoDataFrame(roads, crs=UTM); segs = gpd.GeoDataFrame(segs, crs=UTM); segs["seg_id"] = np.arange(len(segs))
    towns = gpd.GeoDataFrame(towns, crs=UTM)
    st = r[(r.bridge | r.tunnel) & r.highway.isin({c for rd in CFG["road"] for c in rd["classes"]})]
    st = st[st.intersects(roads.buffer(20).union_all())][["highway", "ref", "name_en", "bridge", "tunnel", "geometry"]]
    dams = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="dams").to_crs(UTM)
    o = R / "data/processed/corridor.gpkg"
    if o.exists():
        arch = R / "data/processed/_archive_2026-11-16"; arch.mkdir(exist_ok=True)
        if not (arch / o.name).exists():
            o.rename(arch / o.name)
        else:
            o.unlink()
    roads.to_file(o, layer="roads"); segs.to_file(o, layer="segments"); st.to_file(o, layer="structures")
    dams[["osm_id", "name", "geometry"]].to_file(o, layer="dams"); towns.to_file(o, layer="towns")
    pd.DataFrame(summ).to_csv(R / "outputs/roads_summary.csv", index=False)
    print("segments", len(segs), "structures", len(st), "tunnels:", sorted(st[st.tunnel].name_en.dropna().unique().tolist()))


if __name__ == "__main__":
    main()
