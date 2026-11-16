"""Elevation for the corridor: Copernicus DEM GLO-30 (DSM, TanDEM-X 2011-2015, 1 arc-second) tiles from the open AWS
bucket, mosaicked and clipped to the corridor buffer, in UTM 43N at 30 m.

    python src/02_dem.py

Tiles: https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N{lat}_00_E{lon}_00_DEM/... (1 x 1 degree COGs).
Downloads via .part files and skips tiles already present (reuse). Ocean or missing tiles are skipped.
Writes data/raw/dem/glo30/*.tif, data/processed/dem_glo30_utm.tif (corridor buffer, bilinear, 30 m) and
data/processed/corridor_buffer.gpkg. Other DEMs (FABDEM DTM, AW3D30, NASADEM, HMA 8 m, Sentinel-1 InSAR DEM) are
compared in src/03_dem_compare.py.
"""
import sys
import tomllib
import urllib.request
from pathlib import Path

import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.warp import Resampling, calculate_default_transform, reproject

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
CFG = tomllib.loads((R / "config/corridor.toml").read_text(encoding="utf-8"))
D = R / "data/raw/dem/glo30"; D.mkdir(parents=True, exist_ok=True)
URL = "https://copernicus-dem-30m.s3.amazonaws.com/{n}/{n}.tif"


def tiles(buf):
    w, s, e, n = buf.to_crs(4326).total_bounds
    for la in range(int(s // 1), int(n // 1) + 1):
        for lo in range(int(w // 1), int(e // 1) + 1):
            yield f"Copernicus_DSM_COG_10_N{la:02d}_00_E{lo:03d}_00_DEM"


def get(name):
    f = D / f"{name}.tif"
    if f.exists():
        return f
    try:
        urllib.request.urlretrieve(URL.format(n=name), f.with_suffix(".part"))
    except Exception as e:
        print("skip", name, str(e)[:60], flush=True); return None
    f.with_suffix(".part").rename(f); print("got", name, flush=True)
    return f


def main():
    roads = gpd.read_file(R / "data/processed/corridor.gpkg", layer="roads")
    buf = gpd.GeoDataFrame(geometry=[roads.buffer(CFG["buffer_km"] * 1000).union_all()], crs=roads.crs)
    buf.to_file(R / "data/processed/corridor_buffer.gpkg")
    fs = [f for f in (get(t) for t in tiles(buf)) if f]
    srcs = [rasterio.open(f) for f in fs]
    arr, tr = merge(srcs); prof = srcs[0].profile
    prof.update(height=arr.shape[1], width=arr.shape[2], transform=tr, compress="deflate", tiled=True, BIGTIFF="IF_SAFER")
    tmp = R / "data/processed/dem_glo30_wgs84.tif"
    with rasterio.open(tmp, "w", **prof) as o:
        o.write(arr)
    with rasterio.open(tmp) as src:
        t, w, h = calculate_default_transform(src.crs, roads.crs, src.width, src.height, *src.bounds, resolution=30)
        p = src.profile.copy(); p.update(crs=roads.crs, transform=t, width=w, height=h, nodata=-9999)
        utm = R / "data/processed/dem_glo30_utm_full.tif"
        with rasterio.open(utm, "w", **p) as o:
            reproject(rasterio.band(src, 1), rasterio.band(o, 1), resampling=Resampling.bilinear, dst_nodata=-9999)
    with rasterio.open(utm) as src:
        a, t = mask(src, buf.geometry, crop=True, nodata=-9999)
        p = src.profile.copy(); p.update(height=a.shape[1], width=a.shape[2], transform=t)
    with rasterio.open(R / "data/processed/dem_glo30_utm.tif", "w", **p) as o:
        o.write(a)
    v = a[a > -9000]
    print("tiles", len(fs), "corridor DEM", a.shape, "elev min/max", float(v.min()), float(v.max()))


if __name__ == "__main__":
    main()
