"""Google GEM-15 (Batchu et al. 2026, experimental, CC-BY 4.0, "created by Google"): 15 m surface (dsm) and bare-earth
terrain (dtm) models from sub-metre optical stereo, for the corridor buffer, from Earth Engine.

    python src/02b_gem15.py

Asset projects/climate-and-sustainability/assets/gem-15_v1 (ImageCollection; bands dsm, dsm_std, dsm_percentile_5/50/95,
dtm). Mosaicked and fetched with computePixels in 1024 x 1024 tiles on a UTM 43N 15 m grid aligned to the GLO-30 grid
(every second cell), keeping only tiles that touch the corridor buffer. Tiles via .part files; existing tiles reused.
Writes data/raw/dem/gem15/tiles/<band>_<row>_<col>.tif and data/processed/dem_gem15_{dsm,dtm,dsm_std}_utm.vrt.
Why it matters here: GLO-30 is a surface model (canopy, buildings) at 30 m; GEM-15's dtm is a terrain model at 15 m,
which is what slope, flow routing and HAND need. Reported offsets against lidar are region means of -4.8 to +1.9 m,
so it is checked against GLO-30 and ICESat-2 in src/04_dem_compare.py before use.
"""
import io
import sys
import time
from pathlib import Path

import ee
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
T = R / "data/raw/dem/gem15/tiles"; T.mkdir(parents=True, exist_ok=True)
ASSET = "projects/climate-and-sustainability/assets/gem-15_v1"
BANDS = ["dtm", "dsm", "dsm_std"]
PX, N = 15, 1024


def main():
    ee.Initialize(project="ee-ahmadammargeo")
    buf = gpd.read_file(R / "data/processed/corridor_buffer.gpkg").to_crs(32643)
    geom = buf.geometry.iloc[0]
    w, s, e, n = geom.bounds
    w, n = np.floor(w / 30) * 30, np.ceil(n / 30) * 30
    img = ee.ImageCollection(ASSET).filterBounds(ee.Geometry.Rectangle(list(buf.to_crs(4326).total_bounds))).mosaic().select(BANDS).toFloat()
    span = PX * N
    rows, cols = int(np.ceil((n - s) / span)), int(np.ceil((e - w) / span))
    todo = [(r, c) for r in range(rows) for c in range(cols)
            if geom.intersects(box(w + c * span, n - (r + 1) * span, w + (c + 1) * span, n - r * span))]
    print("tiles", len(todo), "of", rows * cols, flush=True)
    for k, (r, c) in enumerate(todo):
        x0, y0 = w + c * span, n - r * span
        out = T / f"tile_E{int(x0)}_N{int(y0)}.tif"     # named by map position, so a changed grid never reuses a wrong tile
        if out.exists():
            continue
        req = {"expression": img, "fileFormat": "GEO_TIFF", "bandIds": BANDS,
               "grid": {"dimensions": {"width": N, "height": N},
                        "affineTransform": {"scaleX": PX, "shearX": 0, "translateX": x0, "shearY": 0, "scaleY": -PX, "translateY": y0},
                        "crsCode": "EPSG:32643"}}
        for t in range(6):
            try:
                data = ee.data.computePixels(req); break
            except Exception as ex:
                if t == 5:
                    print("failed", r, c, str(ex)[:80], flush=True); data = None
                else:
                    time.sleep(10 * (t + 1))
        if data is None:
            continue
        part = out.with_suffix(".part"); part.write_bytes(data); part.rename(out)
        if k % 20 == 0:
            print(k, "/", len(todo), flush=True)
    tifs = sorted(T.glob("tile_E*_N*.tif"))
    from osgeo import gdal
    for i, b in enumerate(BANDS, 1):
        gdal.BuildVRT(str(R / f"data/processed/dem_gem15_{b}_utm.vrt"), [str(f) for f in tifs], bandList=[i], srcNodata=float("-inf"), VRTNodata=-9999)
    print("done", len(tifs), "tiles")


if __name__ == "__main__":
    main()
