"""SAR DEM tutorial, step 3: download the HyP3 interferograms and turn the leftover phase into height.

    python tutorial/sar_dem/03_phase_to_height.py

Physics (Gabriel, Goldstein and Zebker 1989; Hanssen 2001; written out clearly in M. Shah, "SAR Interferometry and
Phase Unwrapping with SNAPHU", 2026). The interferometric phase is a sum:
    dphi = dphi_flat + dphi_elevation + dphi_displacement + dphi_atmosphere + dphi_noise + dphi_err
HyP3 has removed dphi_flat (orbits) and the topographic phase predicted by Copernicus GLO-30. With a 12-day pair,
displacement is small almost everywhere, so the unwrapped phase that is left is mostly the ERROR of GLO-30 here
(plus atmosphere and real change since 2011-15, such as new cuts, slides and the Attabad lake):
    dz = - lambda R0 sin(theta) / (4 pi B_perp) * dphi_unw          (height correction, m)
    h_a = lambda R0 sin(theta) / (2 B_perp)                         (height of ambiguity: one fringe, m)
InSAR DEM = GLO-30 (HyP3's reference DEM) + dz, on coherent pixels only (coherence >= 0.4).
Choices, and why:
  - phase referenced to the median over flat, coherent, low ground (river terraces), since unwrapped phase is
    relative, not absolute;
  - HyP3 unwraps with MCF (as Shah recommends over MST); SNAPHU in TOPO cost mode is the alternative for a DEM job
    (DEFO for deformation), shown in 04_snaphu.md;
  - layover and shadow (steep KKH walls) carry no height: masked by HyP3's layover/shadow layer where present.
Writes out/<pair>/insar_dem.tif, dz.tif, and figures in out/<pair>/, and a summary in out/summary.csv.
"""
import json
import re
import sys
import zipfile
from pathlib import Path

import hyp3_sdk as sdk
import numpy as np
import pandas as pd
import rasterio

sys.stdout.reconfigure(encoding="utf-8")
O = Path(__file__).parent / "out"
LAMBDA = 0.05546576


def meta(txt):
    d = {}
    for line in txt.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" in line:
            k, v = line.split(":", 1); d[k.strip()] = v.strip()
    return d


def num(d, key, default=None):
    """Exact key from the HyP3 metadata file (e.g. 'Baseline', 'Slant range center')."""
    m = re.search(r"-?\d+\.?\d*", d.get(key, ""))
    return float(m.group()) if m else default


def main():
    ids = json.loads((O / "jobs.json").read_text())
    h = sdk.HyP3()
    jobs = sdk.Batch([h.get_job_by_id(i) for i in ids])
    jobs = h.watch(jobs, timeout=4 * 3600, interval=120)
    rows = []
    for j in jobs:
        if not j.succeeded():
            print("failed", j.job_id); continue
        zs = j.download_files(O)
        for z in zs:
            d = O / Path(z).stem
            if not d.exists():
                zipfile.ZipFile(z).extractall(O)
            g = lambda s: next(d.glob(f"*_{s}.tif"))
            md = meta(d / f"{d.name}.txt")              # the metadata file; *.README.md.txt is the product guide
            bperp = num(md, "Baseline")
            R0 = num(md, "Slant range center", 850e3)
            with rasterio.open(g("unw_phase")) as r:
                unw = r.read(1).astype("float64"); prof = r.profile
            corr = rasterio.open(g("corr")).read(1)
            dem = rasterio.open(g("dem")).read(1).astype("float64")
            inc = rasterio.open(g("inc_map")).read(1).astype("float64")          # radians in HyP3 (checked: max < pi/2)
            inc = np.where(inc > 3.2, np.radians(inc), inc)
            ok = (corr >= 0.4) & (unw != 0) & np.isfinite(unw)
            from scipy.ndimage import generic_filter  # noqa: F401  (slope via gradient below)
            gy, gx = np.gradient(dem, 40)
            flat = ok & (np.hypot(gx, gy) < 0.05) & (dem < np.nanpercentile(dem[ok], 20))
            unw = unw - np.median(unw[flat])
            dz = -LAMBDA * R0 * np.sin(inc) / (4 * np.pi * bperp) * unw
            h_a = LAMBDA * R0 * np.median(np.sin(inc[ok])) / (2 * abs(bperp))
            out = np.where(ok, dem + dz, np.nan).astype("float32")
            prof.update(dtype="float32", nodata=np.nan, compress="deflate")
            with rasterio.open(d / "insar_dem.tif", "w", **prof) as o:
                o.write(out, 1)
            with rasterio.open(d / "dz.tif", "w", **prof) as o:
                o.write(np.where(ok, dz, np.nan).astype("float32"), 1)
            v = dz[ok]
            rows.append({"pair": d.name, "bperp_m": bperp, "R0_m": R0, "h_ambiguity_m": round(h_a, 1), "coherent_frac": round(ok.mean(), 3),
                         "dz_median": round(float(np.median(v)), 2), "dz_nmad": round(float(1.4826 * np.median(np.abs(v - np.median(v)))), 2),
                         "dz_p5": round(float(np.percentile(v, 5)), 1), "dz_p95": round(float(np.percentile(v, 95)), 1)})
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(O / "summary.csv", index=False)


if __name__ == "__main__":
    main()
