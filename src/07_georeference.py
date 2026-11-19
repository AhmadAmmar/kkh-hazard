"""Georeference a scanned map from ground control points, then digitise its colour-coded units into polygons.

    python src/07_georeference.py georef <scan.tif|png|jpg> <gcps.csv> [--order 1|2|tps] [--crs 32643]
    python src/07_georeference.py digitise <georeferenced.tif> <legend.csv> [--min-area-ha 5]

Step 1, georef. gcps.csv has columns px, py (pixel column and row in the scan, from any image viewer or QGIS), lon,
lat (WGS 84, from a recognisable point: a river confluence, a bridge, a peak, a graticule crossing on the map). Use
>= 6 points for a first-order (affine) fit, >= 10 for second order, >= 20 well spread for thin-plate spline (tps,
which bends the scan through every point: only for old, distorted paper maps). The fit is checked by leave-one-out:
each point is predicted from the others; the RMS of those errors (m) is the honest accuracy and is written to
<scan>_georef_report.csv. A map at 1:650,000 drawn to 0.5 mm is at best ~ 325 m accurate; do not expect better.
Step 2, digitise. legend.csv has unit, r, g, b (a colour sampled from each legend box). Every pixel is assigned to
the nearest legend colour in CIELAB space (perceptual distance), unless it is farther than a threshold (text, lines,
hatching: left unassigned), a 5 x 5 majority filter removes speckle, and the units are polygonised and dissolved;
polygons below --min-area-ha are dropped. Output <scan>_units.gpkg with unit and area, ready to check by eye against
the scan and correct by hand where colours are shared or overprinted.
GDAL command-line tools are used (the GDAL Python bindings in this environment are broken).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
BIN = Path("C:/Users/Owner/anaconda3/envs/eo/Library/bin")


def fit(px, py, X, Y, order):
    """Polynomial fit pixel -> map; returns predict function."""
    def design(u, v):
        """Polynomial design matrix of scan coordinates (u, v) to the given order."""
        cols = [np.ones_like(u), u, v]
        if order >= 2:
            cols += [u * u, u * v, v * v]
        return np.column_stack(cols)
    A = design(px, py)
    cx, *_ = np.linalg.lstsq(A, X, rcond=None); cy, *_ = np.linalg.lstsq(A, Y, rcond=None)
    return lambda u, v: (design(u, v) @ cx, design(u, v) @ cy)


def georef(scan, gcps, order, crs):
    """Georeference a scanned map from ground control points: report each point's
    leave-one-out error, then warp with GDAL (polynomial order or thin-plate spline)."""
    g = pd.read_csv(gcps)
    pts = gpd.GeoSeries(gpd.points_from_xy(g.lon, g.lat), crs=4326).to_crs(crs)
    X, Y = pts.x.values, pts.y.values
    o = 1 if order == "tps" else int(order)
    loo = []
    for i in range(len(g)):
        m = np.arange(len(g)) != i
        f = fit(g.px.values[m], g.py.values[m], X[m], Y[m], o)
        px_, py_ = f(g.px.values[i:i + 1], g.py.values[i:i + 1])
        loo.append(float(np.hypot(px_[0] - X[i], py_[0] - Y[i])))
    g["loo_error_m"] = np.round(loo, 1)
    scan = Path(scan)
    g.to_csv(scan.with_name(scan.stem + "_georef_report.csv"), index=False)
    print(f"{len(g)} GCPs, leave-one-out RMS {np.sqrt(np.mean(np.square(loo))):.0f} m (order {order}); worst:",
          g.sort_values("loo_error_m").tail(3)[["px", "py", "loo_error_m"]].values.tolist())
    tmp = scan.with_name(scan.stem + "_gcp.vrt")
    args = [str(BIN / "gdal_translate.exe"), "-of", "VRT", "-a_srs", f"EPSG:{crs}"]
    for (_, r), x, y in zip(g.iterrows(), X, Y):
        args += ["-gcp", str(r.px), str(r.py), str(x), str(y)]
    subprocess.run(args + [str(scan), str(tmp)], check=True)
    out = scan.with_name(scan.stem + "_georef.tif")
    warp = [str(BIN / "gdalwarp.exe"), "-overwrite", "-r", "near", "-co", "COMPRESS=DEFLATE", "-t_srs", f"EPSG:{crs}"]
    warp += ["-tps"] if order == "tps" else ["-order", str(o)]
    subprocess.run(warp + [str(tmp), str(out)], check=True)
    print("written", out)


def lab(rgb):
    """RGB colours converted to CIELAB."""
    from skimage.color import rgb2lab
    return rgb2lab(np.asarray(rgb, dtype="float64").reshape(-1, 1, 3) / 255).reshape(-1, 3)


def digitise(tif, legend, min_area_ha, max_de=18):
    """Polygons of each legend colour in a georeferenced map: pixels within max_de (CIELAB distance) of
        the legend colour, smoothed, polygons below min_area_ha dropped.
    """
    import rasterio
    from rasterio.features import shapes
    from scipy.ndimage import generic_filter
    from shapely.geometry import shape
    leg = pd.read_csv(legend)
    with rasterio.open(tif) as r:
        a = r.read([1, 2, 3]).transpose(1, 2, 0); tr = r.transform; crs = r.crs
    L = lab(a.reshape(-1, 3)); C = lab(leg[["r", "g", "b"]].values)
    d = np.linalg.norm(L[:, None, :] - C[None, :, :], axis=2)
    k = d.argmin(1).astype("int16"); k[d.min(1) > max_de] = -1
    k = k.reshape(a.shape[:2])
    maj = generic_filter(k, lambda v: np.bincount(v[v >= 0].astype(int)).argmax() if (v >= 0).any() else -1, size=5).astype("int16")
    polys = [{"unit": leg.unit.iloc[int(v)], "geometry": shape(gm)} for gm, v in shapes(maj, mask=maj >= 0, transform=tr)]
    g = gpd.GeoDataFrame(polys, crs=crs).dissolve("unit").explode(index_parts=False).reset_index()
    g["area_ha"] = g.area / 1e4
    g = g[g.area_ha >= min_area_ha]
    out = Path(tif).with_name(Path(tif).stem + "_units.gpkg"); g.to_file(out)
    print(g.groupby("unit").area_ha.sum().round(0).sort_values(ascending=False).to_string()); print("written", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--order", default="1"); ap.add_argument("--crs", type=int, default=32643)
    ap.add_argument("--min-area-ha", type=float, default=5)
    x = ap.parse_args()
    georef(x.a, x.b, x.order, x.crs) if x.cmd == "georef" else digitise(x.a, x.b, x.min_area_ha)
