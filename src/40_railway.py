"""A planning-level railway alignment for a mountain line, first case Havelian - Khunjerab (the CPEC proposal), on open
terrain data. Not an engineering design: a reproducible first answer to "where could the line go, how steep, how much
of it in tunnel or on bridge, and roughly what would it cost", to compare options and to argue from.

    python src/40_railway.py [--line havelian_khunjerab] [--grade 2.5] [--cell 90]

Method
1. Horizontal alignment. Least-cost path on a 90 m grid (Copernicus GLO-30 resampled) inside a 15 km band around the
   existing road (the valleys are the only routes), cost = 1 + (slope / 12 deg)^2 + hazard (when the hazard layers
   exist), from the start to the end point; then smoothed so that the curve radius is >= R_min (500 m), by a moving
   average of growing width until the minimum radius holds.
2. Vertical alignment. Dynamic programming along the chainage (100 m steps, 1 m height levels): at each step the rail
   may rise or fall by at most the ruling gradient (default 2.5 %, a mountain main line; 1.2 % is a heavy-freight
   ideal); each step costs open line if the rail is within +/- 12 m of the ground (cut or embankment, cost growing with
   depth), tunnel if the rail is more than 12 m below the ground, bridge or viaduct if more than 12 m above. Rail
   may not be below the maximum pool of reservoirs it passes (Diamer-Bhasha 1160 m, Dasu 950 m, + 15 m freeboard).
   The path of lowest total cost is the profile; it shows directly where height must be won by tunnels and spirals.
3. Report: length, highest point, steepest sustained grade, share and list of tunnels and bridges (> 200 m), and
   indicative cost from config/railway.toml unit costs (ranges, with their sources stated as assumptions).
Writes data/processed/railway_<line>.gpkg (line, structures), outputs/railway_<line>_profile.csv, _summary.csv,
figures/railway_<line>_{map,profile}.png.
"""
import argparse
import sys
import tomllib
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from shapely.geometry import LineString, Point

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
P = R / "data/processed"
CFG = tomllib.loads((R / "config/railway.toml").read_text(encoding="utf-8"))


def dem_grid(cell, bounds):
    with rasterio.open(P / "dem_glo30_utm_full.tif") as r:
        win = rasterio.windows.from_bounds(*bounds, r.transform)
        h, w = int(win.height * r.res[1] / cell), int(win.width * r.res[0] / cell)
        z = r.read(1, window=win, out_shape=(h, w), resampling=Resampling.average).astype("float64")
        tr = rasterio.windows.transform(win, r.transform) * rasterio.Affine.scale(r.width and win.width / w, win.height / h)
    z[z < -100] = np.nan
    return z, tr


def horizontal(line, cell):
    from skimage.graph import route_through_array
    roads = gpd.read_file(P / "corridor.gpkg", layer="roads")
    band = roads[roads.road.isin(line["follow"])].buffer(line.get("band_km", 15) * 1000).union_all()
    a, b = gpd.GeoSeries(gpd.points_from_xy([line["start"][1], line["end"][1]], [line["start"][2], line["end"][2]]), crs=4326).to_crs(roads.crs)
    x0, y0, x1, y1 = band.bounds
    z, tr = dem_grid(cell, (x0, y0, x1, y1))
    gy, gx = np.gradient(np.nan_to_num(z, nan=np.nanmean(z)), cell)
    slope = np.degrees(np.arctan(np.hypot(gx, gy)))
    from rasterio.features import geometry_mask
    inside = ~geometry_mask([band], z.shape, tr)
    cost = 1 + (slope / 12) ** 2
    cost[~inside | ~np.isfinite(z)] = 1e6
    # future reservoirs: ground below the full pool (+ freeboard) near the dam will be lake; the line must run above it
    xs = tr.c + tr.a * (np.arange(z.shape[1]) + 0.5); ys = tr.f + tr.e * (np.arange(z.shape[0]) + 0.5)
    X, Y = np.meshgrid(xs, ys)
    for res in CFG.get("reservoir", []):
        d0 = gpd.GeoSeries(gpd.points_from_xy([res["lon"]], [res["lat"]]), crs=4326).to_crs(roads.crs).iloc[0]
        lake = (np.hypot(X - d0.x, Y - d0.y) < res["reach_km"] * 1000) & (z < res["pool_m"] + res["freeboard_m"])
        cost[lake] = 60.0
    # stay near the valley: penalise climbing far above the road's own height at the nearest road cell (a railway
    # follows valleys on benches along their sides, it does not cross ridges)
    from rasterio.features import rasterize
    from scipy.ndimage import distance_transform_edt
    rr = rasterize([(gm, 1) for gm in roads[roads.road.isin(line["follow"])].geometry], out_shape=z.shape, transform=tr) > 0
    _, (iy, ix) = distance_transform_edt(~rr, return_indices=True)
    zref = np.nan_to_num(z, nan=0)[iy, ix]
    above = np.maximum(0, np.nan_to_num(z, nan=9999) - zref - 50)
    cost += (above / 60) ** 2
    inv = ~tr
    rc = lambda p: tuple(int(v) for v in reversed(inv * (p.x, p.y)))
    path, _ = route_through_array(cost, rc(a), rc(b), fully_connected=True, geometric=True)
    xy = np.array([tr * (c + 0.5, r_ + 0.5) for r_, c in path])
    return xy, roads.crs


def smooth(xy, rmin):
    """Resample to 50 m, then moving-average with growing window until the minimum radius (3-point) holds."""
    ln = LineString(xy); n = int(ln.length // 50)
    p = np.array([ln.interpolate(i * 50).coords[0] for i in range(n + 1)])
    for w in (1, 3, 5, 7, 9, 13, 17, 21, 31, 41, 61):
        k = np.ones(w) / w
        q = np.column_stack([np.convolve(np.pad(p[:, i], w // 2, mode="edge"), k, "valid") for i in (0, 1)])
        a, b, c = q[:-10], q[5:-5], q[10:]                     # radius over 250 m chords (3 points 250 m apart)
        ab, bc, ca = np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1)
        area = np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])) / 2
        rad = ab * bc * ca / np.maximum(4 * area, 1e-9)
        if np.percentile(rad, 2) >= rmin:
            return q, w
    return q, w


def vertical(ground, step, grade, floors, c):
    """DP over chainage x height. Levels every 1 m between min and max feasible heights."""
    lo, hi = np.nanmin(ground) - 400, np.nanmax(ground) + 100
    H = np.arange(np.floor(lo), np.ceil(hi) + 1, 1.0)
    dmax = grade / 100 * step
    k = int(np.floor(dmax))                                  # max level change per step
    n = len(ground)
    def stepcost(i):
        d = H - ground[i]
        tun = c["tunnel"] * (1 + np.maximum(0, -d - c["tunnel_depth_m"]) / 1500)       # deeper: longer access, more heat
        bri = c["bridge"] * (1 + np.maximum(0, d - c["bridge_height_m"]) / 40)          # taller piers cost more, fast
        cst = np.where(d < -c["tunnel_depth_m"], tun, np.where(d > c["bridge_height_m"], bri,
                       c["open"] * (1 + (np.abs(d) / c["tunnel_depth_m"]) ** 2)))
        cst = cst * step / 1000
        cst[d > c["max_bridge_m"]] = 1e9                                                 # no viaduct above ~150 m
        cst[H < floors[i]] = 1e9
        return cst
    acc = stepcost(0).copy(); acc[np.abs(H - ground[0]) > 5] = 1e9
    back = np.zeros((n, len(H)), "int16")
    for i in range(1, n):
        best = np.full(len(H), np.inf); arg = np.zeros(len(H), "int16")
        for s in range(-k, k + 1):
            sh = np.roll(acc, s)
            if s > 0:
                sh[:s] = np.inf
            elif s < 0:
                sh[s:] = np.inf
            better = sh < best; best[better] = sh[better]; arg[better] = s
        acc = best + stepcost(i); back[i] = arg
    # the line may end at the pass or in a summit tunnel under it (portal up to `summit_tunnel_max_m` lower)
    end = np.where((H <= ground[-1] + 5) & (H >= ground[-1] - c.get("summit_drop_m", 0)), acc, np.inf).argmin()
    prof = np.zeros(n, "int64"); prof[-1] = end
    for i in range(n - 1, 0, -1):
        prof[i - 1] = prof[i] - back[i, prof[i]]
    return H[prof], acc[end]


def main(name, grade, cell):
    line = CFG["line"][name]; c = CFG["unit_cost_musd_per_km"]
    xy, crs = horizontal(line, cell)
    q, w = smooth(xy, line.get("r_min_m", 500))
    ln = LineString(q)
    step = 100.0
    ch = np.arange(0, ln.length, step)
    pts = [ln.interpolate(s) for s in ch]
    with rasterio.open(P / "dem_glo30_utm_full.tif") as r:
        ground = np.array([v[0] for v in r.sample([(p.x, p.y) for p in pts])], float)
    ground[ground < -100] = np.nan
    ground = pd.Series(ground).interpolate(limit_direction="both").values
    floors = np.full(len(ch), -1e9)
    for res in CFG.get("reservoir", []):
        poly = gpd.GeoSeries(gpd.points_from_xy([res["lon"]], [res["lat"]]), crs=4326).to_crs(crs).buffer(res["reach_km"] * 1000).iloc[0]
        inr = np.array([poly.contains(p) for p in pts]) & (ground < res["pool_m"] + res["freeboard_m"])
        floors[inr] = res["pool_m"] + res["freeboard_m"]
    cc = dict(open=np.mean(c["open_line"]), tunnel=np.mean(c["tunnel"]), bridge=np.mean(c["bridge"]),
              tunnel_depth_m=CFG["design"]["tunnel_depth_m"], bridge_height_m=CFG["design"]["bridge_height_m"],
              max_bridge_m=CFG["design"].get("max_bridge_m", 150), summit_drop_m=line.get("summit_drop_m", 0))
    rail, total = vertical(ground, step, grade, floors, cc)
    d = rail - ground
    kind = np.where(d < -cc["tunnel_depth_m"], "tunnel", np.where(d > cc["bridge_height_m"], "bridge", "open"))
    prof = pd.DataFrame({"chainage_km": ch / 1000, "x": [p.x for p in pts], "y": [p.y for p in pts], "ground_m": ground.round(1),
                         "rail_m": rail.round(1), "kind": kind, "grade_pct": np.r_[0, np.diff(rail) / step * 100].round(2)})
    O = R / "outputs"; prof.to_csv(O / f"railway_{name}_profile.csv", index=False)
    # structures
    runs = (prof.kind != prof.kind.shift()).cumsum()
    st = prof.groupby(runs).agg(kind=("kind", "first"), start_km=("chainage_km", "first"), end_km=("chainage_km", "last"),
                                max_depth_or_height_m=("rail_m", lambda v: float(np.max(np.abs(v - prof.loc[v.index, "ground_m"])))))
    st["length_km"] = (st.end_km - st.start_km + step / 1000).round(2)
    st = st[(st.kind != "open") & (st.length_km >= 0.2)]
    lens = {k: prof.kind.eq(k).sum() * step / 1000 for k in ("open", "tunnel", "bridge")}
    cost = {k: [lens[kk] * c[{"open": "open_line", "tunnel": "tunnel", "bridge": "bridge"}[kk]][j] for kk in lens] for j, k in ((0, "low"), (1, "high"))}
    sustained = pd.Series(prof.grade_pct).rolling(50, center=True).mean().abs().max()
    summ = {"line": line["name"], "length_km": round(ln.length / 1000, 1), "ruling_grade_pct": grade, "smoothing_window": int(w),
            "start_m": round(ground[0]), "highest_rail_m": round(rail.max()), "end_m": round(ground[-1]),
            "steepest_5km_avg_grade_pct": round(float(sustained), 2),
            **{f"{k}_km": round(v, 1) for k, v in lens.items()}, "tunnels_n": int((st.kind == "tunnel").sum()), "bridges_n": int((st.kind == "bridge").sum()),
            "longest_tunnel_km": float(st[st.kind == "tunnel"].length_km.max() if (st.kind == "tunnel").any() else 0),
            "cost_low_busd": round(sum(cost["low"]) / 1000, 1), "cost_high_busd": round(sum(cost["high"]) / 1000, 1)}
    pd.DataFrame([summ]).to_csv(O / f"railway_{name}_summary.csv", index=False)
    st.to_csv(O / f"railway_{name}_structures.csv", index=False)
    g = gpd.GeoDataFrame([{"name": line["name"], "grade_pct": grade}], geometry=[ln], crs=crs)
    g.to_file(P / f"railway_{name}.gpkg", layer="line")
    for k, v in summ.items():
        print(f"  {k}: {v}")
    figures(name, prof, st, g, line)


def figures(name, prof, st, g, line):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    F = R / "figures"
    fig, ax = plt.subplots(figsize=(14, 4.5))
    ax.fill_between(prof.chainage_km, prof.ground_m.min() - 200, prof.ground_m, color="#d9cbb0", label="ground (GLO-30)")
    ax.plot(prof.chainage_km, prof.rail_m, color="k", lw=1.2, label="rail")
    for k, col in (("tunnel", "#7b3294"), ("bridge", "#008837")):
        m = prof.kind == k
        ax.scatter(prof.chainage_km[m], prof.rail_m[m], s=2, color=col, label=k, zorder=3)
    ax.set_xlabel("chainage (km)"); ax.set_ylabel("height (m)"); ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"{line['name']}: planning-level profile, ruling grade {CFG['design']['ruling_grade_pct']} %", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(F / f"railway_{name}_profile.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    roads = gpd.read_file(P / "corridor.gpkg", layer="roads")
    fig, ax = plt.subplots(figsize=(9, 11))
    roads.plot(ax=ax, color="0.75", lw=0.6)
    g.plot(ax=ax, color="k", lw=1.6)
    for k, col in (("tunnel", "#7b3294"), ("bridge", "#008837")):
        m = prof.kind == k; ax.scatter(prof.x[m], prof.y[m], s=3, color=col, label=k, zorder=3)
    b = g.total_bounds; ax.set_xlim(b[0] - 20e3, b[2] + 20e3); ax.set_ylim(b[1] - 20e3, b[3] + 20e3)
    ax.set_xticks([]); ax.set_yticks([]); ax.legend(); ax.set_title(f"{line['name']}: planning-level alignment")
    fig.savefig(F / f"railway_{name}_map.png", dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="havelian_khunjerab"); ap.add_argument("--grade", type=float, default=None)
    ap.add_argument("--cell", type=float, default=90)
    a = ap.parse_args()
    main(a.line, a.grade or CFG["design"]["ruling_grade_pct"], a.cell)
