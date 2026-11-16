"""Ground: soil, rock, faults and minerals for the study box, from the best open sources, each saved raw then clipped.

    python src/06_geology_soil.py [--only soil,faults,minerals,emit,aster]

Sources and resolution (finest open first):
  soil      ISRIC SoilGrids 2.0 (250 m; Earth Engine projects/soilgrids-isric): clay, sand, silt, coarse fragments
            (cfvo), bulk density (bdod), organic carbon (soc), pH, at 0-5, 5-15, 15-30, 30-60, 60-100, 100-200 cm;
            OpenLandMap USDA texture class (250 m). Mountain soils are thin and poorly sampled here: SoilGrids is a
            model with wide uncertainty; it is used as a covariate, not as truth.
  faults    HimaTibetMap (Styron, Taylor, Yin et al.; active faults, sutures) from GitHub, besides GEM (src/05).
  minerals  USGS Mineral Resources Data System (MRDS) occurrences, clipped.
  emit      NASA EMIT L2B mineral identification (60 m, 2022-now; imaging spectrometer on the ISS): granule list for the
            box from NASA CMR, and download (Earthdata login in ~/_netrc).
  aster     ASTER L1T (15 m VNIR, 30 m SWIR; SWIR valid to April 2008) mineral indices computed in Earth Engine:
            AlOH (clays, micas: b5+b7 / b6), MgOH and carbonate (b7+b9 / b8), ferric iron (b2/b1), ferrous (b5/b3 + b1/b2)
            and the Ninomiya quartz index from TIR (b11^2 / (b10 x b12)), medians of snow-free, cloud-free scenes.
Not scriptable (manual): USGS geologic map of South Asia (geo8ag, ScienceBase blocks scripts: download in a browser
to data/raw/geology/usgs_geo8ag/), and scanned maps that need georeferencing (Searle and Khan 1996, 1:650,000;
GSP KKH road log 2000; Kohistan 1:250,000), done in src/07_georeference.py from ground control points.
Writes data/raw/geology/*, data/processed/ground/*.tif (UTM 43N) and data/processed/inventories.gpkg layers.
"""
import argparse
import io
import json
import netrc
import os
import sys
import time
import tomllib
import urllib.request
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
CFG = tomllib.loads((R / "config/corridor.toml").read_text(encoding="utf-8"))
W, S, E, N = CFG["box"]
RAW, OUT = R / "data/raw/geology", R / "data/processed/ground"
RAW.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
INV = R / "data/processed/inventories.gpkg"


def fetch(url, f, headers=None):
    f = Path(f); f.parent.mkdir(parents=True, exist_ok=True)
    if f.exists():
        return f
    req = urllib.request.Request(url, headers=headers or {"User-Agent": "research-client"})
    with urllib.request.urlopen(req, timeout=600) as r, open(f.with_suffix(f.suffix + ".part"), "wb") as o:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            o.write(b)
    f.with_suffix(f.suffix + ".part").rename(f)
    return f


def ee_grid(img, name, scale):
    """Fetch an Earth Engine image over the box in UTM 43N tiles with computePixels; mosaic to one GeoTIFF."""
    import ee
    import rasterio
    from rasterio.merge import merge
    box = gpd.GeoSeries.from_xy([W, E], [S, N], crs=4326).to_crs(32643).total_bounds
    x0, y1 = np.floor(box[0] / scale) * scale, np.ceil(box[3] / scale) * scale
    n = 1024; span = n * scale
    tiles = []
    d = OUT / f"_tiles_{name}"; d.mkdir(exist_ok=True)
    for r_ in range(int(np.ceil((y1 - box[1]) / span))):
        for c in range(int(np.ceil((box[2] - x0) / span))):
            f = d / f"{r_}_{c}.tif"
            if not f.exists():
                req = {"expression": img, "fileFormat": "GEO_TIFF",
                       "grid": {"dimensions": {"width": n, "height": n}, "crsCode": "EPSG:32643",
                                "affineTransform": {"scaleX": scale, "shearX": 0, "translateX": x0 + c * span, "shearY": 0,
                                                    "scaleY": -scale, "translateY": y1 - r_ * span}}}
                for k in range(5):
                    try:
                        f.with_suffix(".part").write_bytes(ee.data.computePixels(req)); f.with_suffix(".part").rename(f); break
                    except Exception as e:
                        if k == 4:
                            raise
                        time.sleep(10 * (k + 1))
            tiles.append(f)
    srcs = [rasterio.open(t) for t in tiles]
    arr, tr = merge(srcs)
    prof = srcs[0].profile; prof.update(height=arr.shape[1], width=arr.shape[2], transform=tr, compress="deflate", tiled=True)
    descs = srcs[0].descriptions
    with rasterio.open(OUT / f"{name}.tif", "w", **prof) as o:
        o.write(arr)
        for i, dsc in enumerate(descs, 1):
            if dsc:
                o.set_band_description(i, dsc)
    print(name, arr.shape, flush=True)


def soil():
    import ee
    ee.Initialize(project="ee-ahmadammargeo")
    for v in ("clay", "sand", "silt", "cfvo", "bdod", "soc", "phh2o"):
        img = ee.Image(f"projects/soilgrids-isric/{v}_mean").select([f"{v}_0-5cm_mean", f"{v}_15-30cm_mean", f"{v}_60-100cm_mean"]).toFloat()
        ee_grid(img, f"soilgrids_{v}", 250)
    ee_grid(ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").select(["b0", "b30"]).toFloat(), "texture_usda", 250)


def faults():
    d = RAW / "himatibetmap"
    z = fetch("https://github.com/HimaTibetMap/HimaTibetMap/archive/refs/heads/master.zip", d / "HimaTibetMap-master.zip")
    zipfile.ZipFile(z).extractall(d)
    gj = sorted(d.glob("HimaTibetMap-master/geojson/*.geojson"))
    for f in gj:
        g = gpd.read_file(f).cx[W:E, S:N]
        if len(g):
            g.to_file(INV, layer=f"htm_{f.stem}"); print("HimaTibetMap", f.stem, len(g), flush=True)


def minerals():
    z = fetch("https://mrdata.usgs.gov/mrds/mrds-csv.zip", RAW / "usgs_mrds/mrds-csv.zip")
    with zipfile.ZipFile(z) as zz:
        name = [n for n in zz.namelist() if n.endswith(".csv")][0]
        df = pd.read_csv(zz.open(name), low_memory=False, encoding="latin-1")
    df = df[(df.longitude.between(W, E)) & (df.latitude.between(S, N))]
    g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs=4326)
    g.to_file(INV, layer="mrds_minerals"); print("MRDS occurrences", len(g), g.get("commod1", pd.Series()).value_counts().head(12).to_dict(), flush=True)


def emit():
    q = (f"https://cmr.earthdata.nasa.gov/search/granules.json?short_name=EMITL2BMIN&bounding_box={W},{S},{E},{N}"
         f"&page_size=2000&sort_key=-start_date")
    js = json.loads(urllib.request.urlopen(q, timeout=120).read())
    ent = js["feed"]["entry"]
    rows = [{"id": e["title"], "start": e["time_start"][:10], "cloud": e.get("cloud_cover"),
             "url": next((l["href"] for l in e["links"] if l["href"].endswith("_MIN.nc") or "L2B_MIN_" in l["href"] and l["href"].endswith(".nc")), None),
             "poly": e.get("polygons", [[None]])[0][0]} for e in ent]
    t = pd.DataFrame(rows); t.to_csv(RAW / "emit_granules.csv", index=False)
    print("EMIT L2B mineral granules over the box:", len(t), t.start.min() if len(t) else "", "-", t.start.max() if len(t) else "", flush=True)
    if not len(t):
        return
    n = netrc.netrc(os.path.expanduser("~/_netrc")); u, _, p = n.authenticators("urs.earthdata.nasa.gov")
    import requests
    s = requests.Session(); s.auth = (u, p)
    for url in t.url.dropna().head(60):                    # the newest 60 (each ~10-40 MB)
        f = RAW / "emit" / url.split("/")[-1]
        if f.exists():
            continue
        f.parent.mkdir(exist_ok=True)
        r = s.get(url, timeout=600); r.raise_for_status()
        f.with_suffix(".part").write_bytes(r.content); f.with_suffix(".part").rename(f)
    print("EMIT files", len(list((RAW / "emit").glob("*.nc"))), flush=True)


def aster():
    import ee
    ee.Initialize(project="ee-ahmadammargeo")
    roi = ee.Geometry.Rectangle([W, S, E, N])
    col = (ee.ImageCollection("ASTER/AST_L1T_003").filterBounds(roi).filterDate("2000-03-01", "2008-04-01")
           .filter(ee.Filter.lt("CLOUDCOVER", 10)).filter(ee.Filter.calendarRange(6, 10, "month"))
           .filter(ee.Filter.listContains("ORIGINAL_BANDS_PRESENT", "B09")))
    def idx(i):
        b = lambda n: i.select(n).toFloat().max(1)
        snow = b("B01").divide(b("B04")).gt(2.5)                 # bright VNIR, dark SWIR: snow and ice
        out = ee.Image.cat(b("B05").add(b("B07")).divide(b("B06")).rename("AlOH"),
                           b("B07").add(b("B09")).divide(b("B08")).rename("MgOH_carbonate"),
                           b("B02").divide(b("B01")).rename("ferric_iron"),
                           b("B05").divide(b("B3N")).add(b("B01").divide(b("B02"))).rename("ferrous_iron"),
                           b("B11").pow(2).divide(b("B10").multiply(b("B12"))).rename("quartz_ninomiya"))
        return out.updateMask(snow.Not())
    print("ASTER scenes (snow-free season, SWIR era):", col.size().getInfo(), flush=True)
    ee_grid(col.map(idx).median().toFloat(), "aster_mineral_indices", 30)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--only", default="soil,faults,minerals,emit,aster")
    for k in ap.parse_args().only.split(","):
        try:
            globals()[k]()
        except Exception as e:
            print(k, "FAILED:", str(e)[:200], flush=True)
