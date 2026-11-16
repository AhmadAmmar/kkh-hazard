"""Corridor roads, bridges, tunnels and dams from OpenStreetMap (Geofabrik Pakistan extract, GDAL OSM driver).

    python src/01_roads.py

Reads data/raw/osm_corridor.gpkg (made with ogr2ogr from pakistan-latest.osm.pbf, box in config/corridor.toml).
Picks each corridor road by OSM ref, then by name, from trunk/primary/secondary ways; flags bridges and tunnels;
cuts each road into 100 m segments (the unit of risk later). Writes data/processed/corridor.gpkg layers roads,
segments, structures and outputs/roads_summary.csv. UTM 43N (EPSG:32643) for lengths.
"""
import re
import sys
import tomllib
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString
from shapely.ops import substring

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
CFG = tomllib.loads((R / "config/corridor.toml").read_text(encoding="utf-8"))
UTM = 32643


def tag(s, k):
    m = re.search(rf'"{k}"=>"([^"]*)"', s or "")
    return m.group(1) if m else None


def main():
    r = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="roads").to_crs(UTM)
    ot = r.other_tags.fillna("")
    r["ref"] = ot.map(lambda s: tag(s, "ref"))
    r["bridge"] = ot.str.contains('"bridge"=>"yes"|"bridge"=>"viaduct"')
    r["tunnel"] = ot.str.contains('"tunnel"=>"yes"')
    r["name_en"] = ot.map(lambda s: tag(s, "name:en")).fillna(r["name"])
    main_ = r[r.highway.isin(["trunk", "primary", "secondary", "trunk_link", "primary_link", "tertiary"])]
    out = []
    for rd in CFG["road"]:
        refs = {x.replace("-", "").upper() for x in rd.get("refs", [])}
        by_ref = main_.ref.fillna("").str.replace("-", "").str.upper().str.split(";").map(lambda v: bool(refs & set(v)))
        nm = "|".join(map(re.escape, rd.get("names", [])))
        by_name = main_.name_en.fillna("").str.contains(nm, case=False) if nm else False
        sel = main_[by_ref | by_name]
        if rd["id"] == "astore":                       # Astore road has no ref: secondary/tertiary ways in the Astore valley
            sel = main_[main_.name_en.fillna("").str.contains("Astore", case=False)]
        out.append(sel.assign(road=rd["id"], road_name=rd["name"]))
    roads = pd.concat(out).drop_duplicates(subset="osm_id")
    roads = gpd.GeoDataFrame(roads, crs=UTM)
    roads["km"] = roads.length / 1000
    # 100 m segments
    segs = []
    for _, w in roads.iterrows():
        g = w.geometry
        n = max(1, int(np.ceil(g.length / 100)))
        for i in range(n):
            s = substring(g, i * g.length / n, (i + 1) * g.length / n)
            if isinstance(s, LineString):
                segs.append({"road": w.road, "osm_id": w.osm_id, "bridge": w.bridge, "tunnel": w.tunnel, "geometry": s})
    segs = gpd.GeoDataFrame(segs, crs=UTM)
    segs["seg_id"] = np.arange(len(segs))
    st = roads[roads.bridge | roads.tunnel][["road", "name_en", "bridge", "tunnel", "km", "geometry"]]
    dams = gpd.read_file(R / "data/raw/osm_corridor.gpkg", layer="dams").to_crs(UTM)
    o = R / "data/processed/corridor.gpkg"
    roads[["osm_id", "road", "road_name", "highway", "ref", "name_en", "bridge", "tunnel", "km", "geometry"]].to_file(o, layer="roads")
    segs.to_file(o, layer="segments"); st.to_file(o, layer="structures"); dams[["osm_id", "name", "geometry"]].to_file(o, layer="dams")
    s = roads.groupby("road_name").agg(km=("km", "sum"), ways=("osm_id", "size"),
                                       bridge_km=("km", lambda v: v[roads.loc[v.index, "bridge"]].sum()),
                                       tunnel_km=("km", lambda v: v[roads.loc[v.index, "tunnel"]].sum())).round(2)
    s.to_csv(R / "outputs/roads_summary.csv")
    print(s.to_string()); print("segments", len(segs), "dams", len(dams), dams.name.dropna().unique()[:20].tolist())
    print("tunnels:", sorted(roads[roads.tunnel].name_en.dropna().unique().tolist()))


if __name__ == "__main__":
    main()
