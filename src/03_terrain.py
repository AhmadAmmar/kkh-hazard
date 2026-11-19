"""Terrain predictors from the corridor DEM (WhiteboxTools), the base of every hazard model.

    python src/03_terrain.py [dem]        (default data/processed/dem_glo30_utm.tif)

Steps, each written to data/processed/terrain/<name>.tif (30 m, UTM 43N):
  filled       breach depressions (least cost), so water can flow out of pits that are DEM artefacts
  slope        degrees;  aspect  degrees from north;  plan_curv, prof_curv  (1/m)
  relief_1km   max - min elevation within 1 km (local relief, a landslide and rockfall proxy)
  rough        terrain ruggedness index
  d8_acc       D8 flow accumulation (cells);  dinf_sca  D-infinity specific catchment area (m)
  twi          topographic wetness index ln(a / tan b)
  spi          stream power index a tan b (erosion and debris flow proxy)
  streams      channels where accumulation >= 1 km2 (1,111 cells)
  hand         height above nearest drainage (m): flood and debris-flow runout proxy
  dist_stream  horizontal distance to the nearest channel (m)
Note on surface models: GLO-30 is a surface model (DSM); in forest it sits on the canopy, which biases slope and HAND.
src/04_dem_compare.py repeats the key layers on FABDEM (a terrain model) to measure that bias.
"""
import sys
from pathlib import Path

import whitebox

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
T = R / "data/processed/terrain"; T.mkdir(parents=True, exist_ok=True)


def relief(dem, out, size=67):
    """Local relief: max minus min elevation in a 67 x 67 cell (2 km) window, i.e. within ~1 km of each cell."""
    import numpy as np
    import rasterio
    from scipy.ndimage import maximum_filter, minimum_filter
    with rasterio.open(dem) as r:
        a = r.read(1).astype("float32"); prof = r.profile; nod = r.nodata
    m = a == nod
    hi = maximum_filter(np.where(m, -1e4, a), size); lo = minimum_filter(np.where(m, 1e4, a), size)
    out_a = np.where(m, -9999, hi - lo).astype("float32")
    prof.update(dtype="float32", nodata=-9999, compress="deflate")
    with rasterio.open(out, "w", **prof) as o:
        o.write(out_a, 1)
    return 0


def main(dem):
    """Terrain derivatives of the DEM with WhiteboxTools breached DEM, slope, aspect, plan and profile
    curvature, 1 km relief, ruggedness, D8 and D-infinity flow, wetness and stream power indices, streams, height
    above stream and distance to stream."""
    wbt = whitebox.WhiteboxTools(); wbt.set_verbose_mode(False); wbt.set_working_dir(str(T))
    p = lambda n: str(T / f"{n}.tif")
    dem = str(Path(dem).resolve())
    steps = [
        ("filled", lambda: wbt.breach_depressions_least_cost(dem, p("filled"), dist=100, fill=True)),
        ("slope", lambda: wbt.slope(dem, p("slope"), units="degrees")),
        ("aspect", lambda: wbt.aspect(dem, p("aspect"))),
        ("plan_curv", lambda: wbt.plan_curvature(dem, p("plan_curv"))),
        ("prof_curv", lambda: wbt.profile_curvature(dem, p("prof_curv"))),
        ("relief_1km", lambda: relief(dem, p("relief_1km"))),
        ("rough", lambda: wbt.ruggedness_index(dem, p("rough"))),
        ("d8_acc", lambda: wbt.d8_flow_accumulation(p("filled"), p("d8_acc"), out_type="cells")),
        ("dinf_sca", lambda: wbt.d_inf_flow_accumulation(p("filled"), p("dinf_sca"), out_type="sca")),
        ("twi", lambda: wbt.wetness_index(p("dinf_sca"), p("slope"), p("twi"))),
        ("spi", lambda: wbt.stream_power_index(p("dinf_sca"), p("slope"), p("spi"))),
        ("streams", lambda: wbt.extract_streams(p("d8_acc"), p("streams"), threshold=1111)),
        ("hand", lambda: wbt.elevation_above_stream(p("filled"), p("streams"), p("hand"))),
        ("dist_stream", lambda: wbt.euclidean_distance(p("streams"), p("dist_stream"))),
    ]
    for name, f in steps:
        if (T / f"{name}.tif").exists():
            print("have", name, flush=True); continue
        rc = f(); print(name, "ok" if rc == 0 else f"failed {rc}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else R / "data/processed/dem_glo30_utm.tif")
