"""Hazard inventories for the corridor box: faults, earthquakes, landslides, glacial lakes, glaciers.

    python src/05_inventories.py

Sources (open; each saved raw under data/raw/<name>/ and clipped to the box in data/processed/inventories.gpkg):
  faults       GEM Global Active Faults (Styron and Pagani 2020), harmonised GeoJSON, GitHub
  quakes       USGS ComCat, M >= 3, 1900-today, FDSN event service (CSV)
  landslides   NASA Global Landslide Catalog (Kirschbaum et al.), rainfall-triggered, 2007-2020 (CSV export)
  glaciers     RGI 7.0 region 14 South Asia West (NSIDC; needs Earthdata) - fetched in src/06 if login allows
  lakes        glacial lake inventory: HMA glacial lakes (Chen et al. 2021, Zenodo) - fetched in src/06
A source that fails is reported and skipped, not faked. Downloads go through .part files and are reused.
"""
import io
import sys
import tomllib
import urllib.request
from datetime import date
from pathlib import Path

import geopandas as gpd
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
CFG = tomllib.loads((R / "config/corridor.toml").read_text(encoding="utf-8"))
W, S, E, N = CFG["box"]
OUT = R / "data/processed/inventories.gpkg"


def fetch(url, f):
    """Download a URL to a file unless it exists; returns the path."""
    f = Path(f); f.parent.mkdir(parents=True, exist_ok=True)
    if f.exists():
        return f
    req = urllib.request.Request(url, headers={"User-Agent": "research-client"})
    with urllib.request.urlopen(req, timeout=300) as r, open(f.with_suffix(".part"), "wb") as o:
        o.write(r.read())
    f.with_suffix(".part").rename(f)
    return f


def save(gdf, layer):
    """Write a layer clipped to the study box to the inventory file."""
    gdf = gdf.cx[W:E, S:N]
    gdf.to_file(OUT, layer=layer)
    print(layer, len(gdf), flush=True)
    return gdf


def main():
    """Fetch the hazard inventories (GEM active faults, USGS earthquakes M3+ since 1900, NASA
    landslide catalogue and others) and save them clipped to the corridor."""
    raw = R / "data/raw"
    try:
        f = fetch("https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/master/geojson/gem_active_faults_harmonized.geojson",
                  raw / "gem_faults/gem_active_faults_harmonized.geojson")
        save(gpd.read_file(f), "faults")
    except Exception as e:
        print("faults failed", str(e)[:100])
    try:
        url = (f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1900-01-01&endtime={date.today()}"
               f"&minlatitude={S - 1}&maxlatitude={N + 1}&minlongitude={W - 1}&maxlongitude={E + 1}&minmagnitude=3")
        f = fetch(url, raw / f"usgs_quakes/quakes_{date.today()}.csv")
        q = pd.read_csv(f)
        save(gpd.GeoDataFrame(q, geometry=gpd.points_from_xy(q.longitude, q.latitude), crs=4326), "quakes")
    except Exception as e:
        print("quakes failed", str(e)[:100])
    try:
        f = fetch("https://data.nasa.gov/api/views/dd9e-wu2v/rows.csv?accessType=DOWNLOAD", raw / "nasa_glc/glc.csv")
        g = pd.read_csv(f)
        lon = [c for c in g.columns if c.lower() in ("longitude", "lon")][0]; lat = [c for c in g.columns if c.lower() in ("latitude", "lat")][0]
        save(gpd.GeoDataFrame(g, geometry=gpd.points_from_xy(g[lon], g[lat]), crs=4326), "landslides_nasa")
    except Exception as e:
        print("NASA landslides failed", str(e)[:100])


if __name__ == "__main__":
    main()
