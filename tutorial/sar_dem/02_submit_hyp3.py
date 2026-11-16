"""SAR DEM tutorial, step 2: ask ASF HyP3 to make the interferograms (GAMMA software, on ASF's servers).

    python tutorial/sar_dem/02_submit_hyp3.py [n_pairs]      (default 2: the best ascending and descending pairs)

Uses the Earthdata login in ~/_netrc (never written here). Each job: INSAR_GAMMA, 10 x 2 looks (~40 m pixels),
with the unwrapped phase, coherence, the look vectors, the incidence angle and the reference DEM (Copernicus GLO-30)
included, because step 3 needs them to turn the leftover phase into height. HyP3 removes the topographic phase
predicted by GLO-30, so the unwrapped phase that remains is: DEM error x (4 pi B_perp / (lambda R sin theta)) +
ground motion + atmosphere + noise. Over 12 days motion is small almost everywhere except active slopes; that is
both the noise in our DEM and the landslide signal we want elsewhere. Job ids are saved to out/jobs.json.
"""
import json
import sys
from pathlib import Path

import hyp3_sdk as sdk
import pandas as pd

O = Path(__file__).parent / "out"


def main(n=2):
    p = pd.read_csv(O / "pairs.csv")
    pick = pd.concat([p[p["dir"] == d].head(n // 2 or 1) for d in ("ASCENDING", "DESCENDING")])
    h = sdk.HyP3()                                    # reads ~/_netrc
    jobs = sdk.Batch()
    for _, r in pick.iterrows():
        jobs += h.submit_insar_job(r.ref, r.sec, name="kkh_dem", looks="10x2", include_dem=True, include_look_vectors=True,
                                   include_inc_map=True, include_los_displacement=False, apply_water_mask=False,
                                   include_displacement_maps=False, phase_filter_parameter=0.6)
        print("submitted", r.d1, r.d2, r["dir"], "bperp", r.bperp_m, flush=True)
    (O / "jobs.json").write_text(json.dumps([j.job_id for j in jobs], indent=1))
    print(h.check_credits() if hasattr(h, "check_credits") else "")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
