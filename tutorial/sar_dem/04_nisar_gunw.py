"""SAR DEM tutorial, step 4: NISAR L-band interferograms (GUNW) over the corridor: why they are poor for heights and
good for motion, and how their coherence compares with Sentinel-1.

    python tutorial/sar_dem/04_nisar_gunw.py

Reads data/raw/nisar/*GUNW*.h5 (downloaded from ASF with the Earthdata login). Per product:
  - perpendicular baseline from the metadata cube (metadata/radarGrid/perpendicularBaseline, by height layer);
  - height of ambiguity h_a = lambda R sin(theta) / (2 |B_perp|): NISAR keeps its orbit in a tight tube, so B_perp is
    a few metres to a few tens of metres and h_a is kilometres: almost no topographic signal;
  - coherence (unwrappedInterferogram/HH/coherenceMagnitude, 80 m) inside the 10 km corridor buffer: share >= 0.4;
  - line-of-sight displacement d = -lambda / (4 pi) x unwrapped phase (sign as in the NISAR product spec: positive
    phase = range increase, i.e. motion away), minus the ionosphere screen, referenced to the corridor median, only in
    the largest connected component and where coherence >= 0.4;
  - sampled at the 100 m road segments: coherence and |d| per road.
Writes tutorial/sar_dem/out/nisar_summary.csv, out/nisar_segments.gpkg, and figures/maps/nisar_*.png.
"""
import glob
import sys
from pathlib import Path

import geopandas as gpd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import from_origin

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[2]
O = Path(__file__).parent / "out"
FIG = R / "figures/maps"
C = 299_792_458.0


def main():
    """Find NISAR interferogram products (GUNW) over the corridor segments."""
    buf = gpd.read_file(R / "data/processed/corridor_buffer.gpkg").to_crs(32643)
    segs = gpd.read_file(R / "data/processed/corridor.gpkg", layer="segments").to_crs(32643)
    mid = segs.geometry.interpolate(0.5, normalized=True)
    rows, sample = [], {}
    for fn in sorted(glob.glob(str(R / "data/raw/nisar/*GUNW*.h5"))):
        name = Path(fn).stem
        tag = "asc" if "_A_" in name else "desc"
        with h5py.File(fn, "r") as f:
            fa = f["science/LSAR/GUNW/grids/frequencyA"]
            lam = C / fa["centerFrequency"][()]
            g = fa["unwrappedInterferogram"]
            x, y = g["xCoordinates"][:], g["yCoordinates"][:]
            unw = g["HH/unwrappedPhase"][:].astype("float64")
            coh = g["HH/coherenceMagnitude"][:].astype("float64")
            cc = g["HH/connectedComponents"][:]
            iono = g["HH/ionospherePhaseScreen"][:].astype("float64")
            rg = f["science/LSAR/GUNW/metadata/radarGrid"]
            bp = rg["perpendicularBaseline"][:]; inc = rg["incidenceAngle"][:]; sr = rg["referenceSlantRange"][:]
            hts = rg["heightAboveEllipsoid"][:]
        k = int(np.argmin(np.abs(hts - 3000)))                     # layer near corridor heights
        B, th, Rr = np.nanmedian(bp[k]), np.radians(np.nanmedian(inc[k])), np.nanmedian(sr[k])
        h_a = lam * Rr * np.sin(th) / (2 * abs(B))
        tr = from_origin(x[0] - 40, y[0] + 40, 80, 80)
        inside = ~geometry_mask(buf.geometry, unw.shape, tr)
        valid = inside & np.isfinite(coh) & (coh > 0)
        big = np.bincount(cc[valid & (cc > 0)]).argmax() if (valid & (cc > 0)).any() else 0
        good = valid & (coh >= 0.4) & (cc == big)
        d = -lam / (4 * np.pi) * (unw - np.where(np.isfinite(iono), iono, 0))
        d = np.where(good, d - np.nanmedian(d[good]), np.nan) * 100      # cm
        rows.append({"product": name[:60], "pass": tag, "lambda_cm": round(lam * 100, 2), "bperp_m": round(float(B), 2),
                     "slant_range_km": round(Rr / 1000, 1), "incidence_deg": round(np.degrees(th), 1), "h_ambiguity_m": round(h_a),
                     "corridor_cells": int(valid.sum()), "coh_median": round(float(np.median(coh[valid])), 3),
                     "share_coh_ge_0.4": round(float((coh[valid] >= 0.4).mean()), 3), "share_used": round(float(good.sum() / valid.sum()), 3),
                     "d_p5_cm": round(float(np.nanpercentile(d, 5)), 1), "d_p95_cm": round(float(np.nanpercentile(d, 95)), 1)})
        print(rows[-1], flush=True)
        col = ((mid.x.values - x[0] + 40) // 80).astype(int); row = ((y[0] + 40 - mid.y.values) // 80).astype(int)
        ok = (row >= 0) & (row < unw.shape[0]) & (col >= 0) & (col < unw.shape[1])
        cs = np.full(len(mid), np.nan); ds = np.full(len(mid), np.nan)
        cs[ok] = coh[row[ok], col[ok]]; ds[ok] = d[row[ok], col[ok]]
        sample[f"coh_{tag}"] = cs; sample[f"dlos_{tag}_cm"] = ds
        # figures
        ext = [x[0], x[-1], y[-1], y[0]]
        w, s_, e, n = buf.total_bounds
        for arr, lab, cm, lim, fname in ((coh, "coherence", "magma", (0, 1), f"nisar_coherence_{tag}"),
                                         (d, "line-of-sight displacement (cm), + = away from satellite", "RdBu_r", (-6, 6), f"nisar_dlos_{tag}")):
            fig, ax = plt.subplots(figsize=(10, 9))
            im = ax.imshow(np.where(inside, arr, np.nan), cmap=cm, vmin=lim[0], vmax=lim[1], extent=ext, interpolation="nearest")
            segs.plot(ax=ax, color="k", lw=0.3)
            ax.set_xlim(w, e); ax.set_ylim(s_, n); ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(im, ax=ax, shrink=0.6, label=lab)
            ax.set_title(f"NISAR L band, {tag}ending, 22 Sep - 4 Oct 2026: {lab.split(' (')[0]}\n"
                         f"B_perp {B:.1f} m, height of ambiguity {h_a / 1000:.1f} km, {100 * rows[-1]['share_coh_ge_0.4']:.0f} % of corridor >= 0.4",
                         fontsize=10)
            fig.text(0.02, 0.02, "NISAR L2 GUNW (provisional), 80 m, ASF. Ionosphere screen removed; largest connected component; "
                     "coherence >= 0.4. UTM 43N. Roads: OSM.", fontsize=7, color="0.4")
            fig.savefig(FIG / f"{fname}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    pd.DataFrame(rows).to_csv(O / "nisar_summary.csv", index=False)
    out = segs[["seg_id", "road", "geometry"]].copy()
    for k_, v in sample.items():
        out[k_] = v
    out.to_file(O / "nisar_segments.gpkg")
    t = out.drop(columns="geometry").groupby("road").agg(**{f"{c}_med": (c, "median") for c in sample if c.startswith("coh")},
                                                         **{f"{c}_share_ge_0.4": (c, lambda v: (v >= 0.4).mean()) for c in sample if c.startswith("coh")})
    print(t.round(3).to_string())
    t.round(3).to_csv(O / "nisar_coherence_by_road.csv")


if __name__ == "__main__":
    main()
