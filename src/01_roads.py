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


def graph(ways):
    G = nx.Graph()
    for idx, g in zip(ways.index, ways.geometry):
        lines = [g] if g.geom_type == "LineString" else list(g.geoms)
        for ln in lines:
            c = [(round(x, 2), round(y, 2)) for x, y in ln.coords]
            for a, b in zip(c, c[1:]):
                if a != b:
                    G.add_edge(a, b, w=float(np.hypot(a[0] - b[0], a[1] - b[1])), way=idx)
    return G


def main():
    r = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="roads").to_crs(UTM)
    ot = r.other_tags.fillna("")
    r["ref"] = ot.map(lambda s: tag(s, "ref"))
    r["bridge"] = ot.str.contains('"bridge"=>"yes"|"bridge"=>"viaduct"')
    r["tunnel"] = ot.str.contains('"tunnel"=>"yes"')
    r["name_en"] = ot.map(lambda s: tag(s, "name:en")).fillna(r["name"])
    roads, segs, towns, summ = [], [], [], []
    for rd in CFG["road"]:
        refs = {x.replace("-", "").upper() for x in rd.get("refs", [])}
        by_ref = r.ref.fillna("").str.replace("-", "").str.upper().str.split(";").map(lambda v: bool(refs & set(v)))
        ways = r[r.highway.isin(rd["classes"]) | by_ref]          # allowed classes, plus any way carrying the road's number
        G = graph(ways)
        comp = max(nx.connected_components(G), key=len)          # keep the main network; islands cannot be routed
        nodes = np.array(list(comp))
        pts = gpd.GeoSeries(gpd.points_from_xy([v[1] for v in rd["via"]], [v[2] for v in rd["via"]]), crs=4326).to_crs(UTM)
        snapped = []
        for (nm, _, _), p in zip(rd["via"], pts):
            d = np.hypot(nodes[:, 0] - p.x, nodes[:, 1] - p.y); i = int(d.argmin())
            snapped.append(tuple(nodes[i]))
            towns.append({"road": rd["id"], "town": nm, "snap_m": round(float(d[i])), "geometry": Point(nodes[i])})
            if d[i] > 1000:
                print(f"  ! {rd['id']}: {nm} is {d[i] / 1000:.1f} km from the {'/'.join(rd['classes'][:2])} network", flush=True)
        used = []
        for k, (a, b) in enumerate(zip(snapped, snapped[1:])):
            try:
                path = nx.shortest_path(G, a, b, weight="w")
            except nx.NetworkXNoPath:
                print(f"  ! {rd['id']}: no path {rd['via'][k][0]} -> {rd['via'][k + 1][0]}", flush=True); continue
            line = LineString(path)
            stretch = f"{rd['via'][k][0]} - {rd['via'][k + 1][0]}"
            ways_used = {G[u][v]["way"] for u, v in zip(path, path[1:])}
            used += list(ways_used)
            br = r.loc[list(ways_used)]
            roads.append({"road": rd["id"], "road_name": rd["name"], "stretch": stretch, "km": line.length / 1000,
                          "bridge_km": br[br.bridge].length.sum() / 1000, "tunnel_km": br[br.tunnel].length.sum() / 1000, "geometry": line})
            n = max(1, int(np.ceil(line.length / 100)))
            sub = br[br.bridge | br.tunnel]
            for i in range(n):
                s = substring(line, i * line.length / n, (i + 1) * line.length / n)
                hit = sub[sub.distance(s.interpolate(0.5, normalized=True)) < 15] if len(sub) else sub
                segs.append({"road": rd["id"], "stretch": stretch, "bridge": bool(len(hit) and hit.bridge.any()),
                             "tunnel": bool(len(hit) and hit.tunnel.any()), "geometry": s})
        rr = [x for x in roads if x["road"] == rd["id"]]
        summ.append({"road": rd["name"], "km": round(sum(x["km"] for x in rr), 1), "stretches": len(rr),
                     "bridge_km": round(sum(x["bridge_km"] for x in rr), 2), "tunnel_km": round(sum(x["tunnel_km"] for x in rr), 2),
                     "worst_snap_m": max(t["snap_m"] for t in towns if t["road"] == rd["id"]),
                     "worst_snap_town": max((t for t in towns if t["road"] == rd["id"]), key=lambda t: t["snap_m"])["town"]})
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
