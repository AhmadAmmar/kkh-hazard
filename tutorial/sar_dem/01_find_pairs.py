"""SAR DEM tutorial, step 1: find Sentinel-1 pairs good for topography over a site on the KKH.

    python tutorial/sar_dem/01_find_pairs.py [lon lat]          (default: Attabad lake, Hunza, 74.82 E 36.32 N)

What makes a pair good for a DEM (not for deformation):
  - a LARGE perpendicular baseline (B_perp): the height of ambiguity h_a = lambda R sin(theta) / (2 B_perp) is the
    height change that wraps the phase once; Sentinel-1's orbital tube is small (B_perp mostly < 150 m), so h_a is
    rarely below ~100 m. That is why Sentinel-1 makes coarse, noisy DEMs and TanDEM-X (B_perp ~ 100-500 m, single pass,
    no time gap) makes good ones. We still do it, to see the physics and to check GLO-30 for change (landslides,
    new cuts, the Attabad lake) since 2011-2015.
  - a SHORT time gap (6 or 12 days) and cold, dry, snow-free conditions: coherence stays high; ice and snow decorrelate.
  - same track and frame (relative orbit), so the geometry matches.
Searches the Alaska Satellite Facility catalogue (no login needed to search), builds the baseline stack for a
reference scene, and ranks pairs by |B_perp| / days. Writes tutorial/sar_dem/out/pairs.csv.
"""
import sys
from pathlib import Path

import asf_search as asf
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
O = Path(__file__).parent / "out"; O.mkdir(exist_ok=True)
LAMBDA = 0.05546576   # m, Sentinel-1 C band
R_SLANT, THETA = 850e3, 0.68   # rough slant range (m) and incidence (rad, ~39 deg) for IW mid-swath


def main(lon=74.82, lat=36.32):
    """Sentinel-1 IW SLC scenes over a point (Oct 2023-Mar 2025) from ASF, paired by track for
        interferometry, written to pairs.csv.
    """
    wkt = f"POINT({lon} {lat})"
    res = asf.geo_search(platform=asf.PLATFORM.SENTINEL1, processingLevel=asf.PRODUCT_TYPE.SLC,
                         beamMode=asf.BEAMMODE.IW, intersectsWith=wkt, start="2023-10-01", end="2025-03-31")
    print("SLC scenes", len(res))
    rows = []
    tracks = pd.Series([r.properties["pathNumber"] for r in res]).value_counts()
    for path in tracks.index[:3]:
        ref = [r for r in res if r.properties["pathNumber"] == path][0]
        stack = ref.stack()                           # all scenes of this track/frame, with baselines to ref
        st = pd.DataFrame([{"scene": s.properties["sceneName"], "date": s.properties["startTime"][:10],
                            "bperp": s.properties.get("perpendicularBaseline"), "dir": s.properties["flightDirection"]}
                           for s in stack]).dropna(subset=["bperp"])
        st["date"] = pd.to_datetime(st.date); st = st.sort_values("date")
        st = st[(st.date >= "2023-10-01") & (st.date <= "2025-03-31")]
        for i in range(len(st)):
            for j in range(i + 1, len(st)):
                dt = (st.date.iloc[j] - st.date.iloc[i]).days
                if dt > 24:
                    break
                bp = st.bperp.iloc[j] - st.bperp.iloc[i]
                rows.append({"path": path, "dir": st.dir.iloc[i], "ref": st.scene.iloc[i], "sec": st.scene.iloc[j],
                             "d1": st.date.iloc[i].date(), "d2": st.date.iloc[j].date(), "days": dt, "bperp_m": round(bp, 1),
                             "h_ambiguity_m": round(LAMBDA * R_SLANT * 0.63 / (2 * abs(bp)), 1) if bp else None,
                             "winter": st.date.iloc[i].month in (11, 12, 1, 2)})
    p = pd.DataFrame(rows)
    p["score"] = p.bperp_m.abs() / p.days
    p = p.sort_values("score", ascending=False)
    p.to_csv(O / "pairs.csv", index=False)
    print(p.head(12).to_string(index=False))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(*(map(float, a))) if a else main()
