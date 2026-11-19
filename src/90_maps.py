"""Maps of every corridor raster made so far, for reading the ground: topography, hydrology, structure, exposure.

    python src/90_maps.py [--area all|attabad|raikot|babusar|skardu|gilgit] [--only name,name]

Each map: the layer over a hillshade (multi-directional, from the same DEM), roads (OSM, by corridor road), tunnels,
dams, active faults (GEM) and M>=4 earthquakes (USGS), scale bar, north arrow, colour bar with units, and a caption
saying what the layer is, how it was made and how to read it. Large rasters are read decimated to the figure size,
so a whole-corridor map and a 15 m zoom come from the same code. Colour maps: perceptual (cividis, viridis, mako-like)
for continuous data, diverging (RdBu) centred on zero for differences, terrain tints for height.
Writes figures/maps/<layer>_<area>.png and a contact sheet figures/maps/index_<area>.png.
"""
import argparse
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LightSource, LinearSegmentedColormap
from rasterio.enums import Resampling
from rasterio.windows import from_bounds

sys.stdout.reconfigure(encoding="utf-8")
R = Path(__file__).resolve().parents[1]
P, T = R / "data/processed", R / "data/processed/terrain"
FIG = R / "figures/maps"; FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
UTM = 32643
# zoom areas: centre lon, lat, half-width km
AREAS = {"all": None, "attabad": (74.82, 36.32, 9), "raikot": (74.60, 35.47, 10), "babusar": (74.04, 35.15, 10),
         "skardu": (75.45, 35.55, 15), "gilgit": (74.33, 35.92, 12), "dasu": (73.25, 35.28, 12)}
ROADC = {"n35": "#111111", "n15": "#7b2d8b", "astore": "#b5651d", "ghizer": "#1f6f8b", "skardu": "#2e7d32"}
TINT = LinearSegmentedColormap.from_list("tint", ["#5a8f5a", "#a7c27a", "#e8dfa8", "#c9a36b", "#9a7556", "#b9b2ac", "#ffffff"])

LAYERS = {   # name: (file, label, cmap, range rule, caption)
    "elevation": (P / "dem_gem15_dtm_utm.vrt", "Height above sea level (m)", TINT, "p1-99",
                  "Google GEM-15 bare-earth terrain model (15 m, from sub-metre stereo; experimental). The Indus, Gilgit and Hunza "
                  "valleys cut 4-6 km below the peaks: the road lives on the valley walls."),
    "elevation_glo30": (P / "dem_glo30_utm_full.tif", "Height above sea level (m)", TINT, "p1-99",
                        "Copernicus GLO-30 surface model (30 m, TanDEM-X 2011-15): includes trees and buildings."),
    "slope": (T / "slope.tif", "Slope (degrees)", "magma_r", (0, 60),
              "Steepness from GLO-30. Above ~35 degrees loose material is near its angle of rest: rockfall and debris slopes."),
    "relief_1km": (T / "relief_1km.tif", "Local relief within ~1 km (m)", "cividis", "p1-99",
                   "Highest minus lowest point in a 2 km window: how much energy a falling rock or debris flow can gain."),
    "hand": (T / "hand.tif", "Height above nearest drainage (m)", "Blues_r", (0, 100),
             "How far above the nearest channel each cell is (channels: catchment >= 1 km2). Low values are where floods and "
             "debris flows run; road on low HAND is exposed."),
    "twi": (T / "twi.tif", "Topographic wetness index", "viridis", "p2-98",
            "ln(upslope area / tan slope): where water gathers. High in valley floors and hollows, the source zones of shallow slides."),
    "spi": (T / "spi.tif", "Stream power index (log10)", "inferno", "log",
            "Upslope area x tan slope: erosive power of flowing water; high where gullies and debris-flow channels cut."),
    "plan_curv": (T / "plan_curv.tif", "Plan curvature", "RdBu_r", "sym",
                  "Across-slope curvature: hollows (convergent, collect water and debris) against noses (divergent)."),
    "dsm_minus_dtm": (None, "Surface minus terrain (m)", "YlGn", (0, 25),
                      "GEM-15 surface minus its bare-earth terrain: tree and building height. Forest belts on the Kaghan and "
                      "Astore sides, orchards in Hunza and Gilgit, towns."),
    "gem15_minus_glo30": (None, "GEM-15 minus GLO-30 (m)", "RdBu_r", (-30, 30),
                          "Difference of the two DEMs: blue where GEM-15 is lower (trees removed, or change since 2011-15 such as new "
                          "cuts and slides), red where higher (new fills, debris, the Attabad lake surface)."),
}


def window_bounds(area, ref):
    """Map bounds of a named area (centre and size in km) in the reference CRS; None for the full corridor."""
    if AREAS[area] is None:
        return None
    lon, lat, km = AREAS[area]
    pt = gpd.GeoSeries(gpd.points_from_xy([lon], [lat]), crs=4326).to_crs(ref).iloc[0]
    return (pt.x - km * 1000, pt.y - km * 1000, pt.x + km * 1000, pt.y + km * 1000)


def read(path, bounds, maxpx=1400):
    """A raster read over the bounds, decimated to at most maxpx pixels, with its extent."""
    with rasterio.open(path) as r:
        win = from_bounds(*bounds, r.transform) if bounds else None
        w = r.width if win is None else int(win.width); h = r.height if win is None else int(win.height)
        f = max(1, max(w, h) / maxpx)
        a = r.read(1, window=win, out_shape=(int(h / f), int(w / f)), resampling=Resampling.average, boundless=bool(win),
                   fill_value=np.nan).astype("float64")
        b = rasterio.windows.bounds(win, r.transform) if win is not None else r.bounds
        nod = r.nodata
    if nod is not None:
        a[a == nod] = np.nan
    a[~np.isfinite(a)] = np.nan; a[a < -500] = np.nan
    return a, tuple(b)


def hillshade(dem, px):
    """Hillshade of a DEM lit from the north-west (half), west and north (a quarter each)."""
    ls = LightSource(azdeg=315, altdeg=40)
    d = np.where(np.isfinite(dem), dem, np.nanmean(dem))
    hs = sum(ls.hillshade(d, vert_exag=1.5, dx=px, dy=px) * w for ls, w in
             ((LightSource(315, 40), 0.5), (LightSource(270, 40), 0.25), (LightSource(0, 40), 0.25)))
    return np.where(np.isfinite(dem), hs, np.nan)


def rng(a, rule):
    """Colour range of an array: a fixed (min, max), 'sym' (+-98th percentile of |a|),
    'log' (2nd-98th percentile) or 'pLO-HI' percentiles."""
    v = a[np.isfinite(a)]
    if isinstance(rule, tuple):
        return rule
    if rule == "sym":
        m = np.percentile(np.abs(v), 98); return (-m, m)
    if rule == "log":                                   # values are already log10 here (taken in main)
        return tuple(np.percentile(v, [2, 98]))
    lo, hi = map(float, rule[1:].split("-"))
    return tuple(np.percentile(v, [lo, hi]))


def overlays(ax, bounds, vec):
    """Draw the vector layers clipped to the bounds."""
    w, s, e, n = bounds
    for name, g in vec.items():
        g = g.cx[w:e, s:n]
        if g.empty:
            continue
        if name == "roads":
            for rd, gg in g.groupby("road"):
                gg.plot(ax=ax, color="white", lw=2.6, zorder=4); gg.plot(ax=ax, color=ROADC.get(rd, "k"), lw=1.3, zorder=5)
        elif name == "tunnels":
            g.plot(ax=ax, color="#ffcc00", lw=3, zorder=6)
        elif name == "faults":
            g.plot(ax=ax, color="#d62728", lw=1.1, ls="--", zorder=6)
        elif name == "quakes":
            q = g[g.mag >= 4]
            ax.scatter(q.geometry.x, q.geometry.y, s=(q.mag - 3.5) ** 2.2 * 6, facecolor="none", edgecolor="#d62728", lw=0.7, zorder=7)
        elif name == "dams":
            ax.scatter(g.geometry.centroid.x, g.geometry.centroid.y, marker="s", s=24, color="#0b4f8a", edgecolor="white", zorder=8)
    ax.set_xlim(w, e); ax.set_ylim(s, n)


def furniture(ax, bounds):
    """Scale bar and north arrow."""
    w, s, e, n = bounds
    span = (e - w) / 1000
    km = [k for k in (1, 2, 5, 10, 20, 50, 100) if k <= span / 4][-1]
    x0, y0 = w + (e - w) * 0.05, s + (n - s) * 0.05
    ax.plot([x0, x0 + km * 1000], [y0, y0], color="k", lw=3, solid_capstyle="butt", zorder=9)
    ax.plot([x0, x0 + km * 1000], [y0, y0], color="white", lw=1, solid_capstyle="butt", zorder=10)
    ax.text(x0 + km * 500, y0 + (n - s) * 0.015, f"{km} km", ha="center", fontsize=8, zorder=10,
            bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
    ax.annotate("N", xy=(0.95, 0.97), xytext=(0.95, 0.88), xycoords="axes fraction", ha="center", fontweight="bold",
                arrowprops=dict(arrowstyle="-|>", color="k"), zorder=10)
    ax.set_xticks([]); ax.set_yticks([])


def legend(fig):
    """Figure legend of the road routes."""
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, lw=2, label=l) for c, l in ((ROADC["n35"], "KKH N-35"), (ROADC["n15"], "Babusar N-15"),
         (ROADC["astore"], "Astore"), (ROADC["ghizer"], "Ghizer"), (ROADC["skardu"], "Skardu S-1"))]
    h += [Line2D([], [], color="#ffcc00", lw=3, label="tunnel"), Line2D([], [], color="#d62728", lw=1, ls="--", label="active fault (GEM)"),
          Line2D([], [], marker="o", mfc="none", mec="#d62728", ls="", label="earthquake M>=4 (USGS)"),
          Line2D([], [], marker="s", color="#0b4f8a", ls="", label="dam (OSM)")]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=7.5, bbox_to_anchor=(0.45, 0.0))


def main(area, only):
    """Hazard maps of an area: hillshaded layers with roads, tunnels, dams and other overlays."""
    vec = {"roads": gpd.read_file(P / "corridor.gpkg", layer="roads").to_crs(UTM),
           "tunnels": gpd.read_file(P / "corridor.gpkg", layer="structures").to_crs(UTM).query("tunnel"),
           "dams": gpd.read_file(P / "corridor.gpkg", layer="dams").to_crs(UTM),
           "faults": gpd.read_file(P / "inventories.gpkg", layer="faults").to_crs(UTM),
           "quakes": gpd.read_file(P / "inventories.gpkg", layer="quakes").to_crs(UTM)}
    if area == "all":
        b0 = tuple(vec["roads"].total_bounds + np.array([-15e3, -15e3, 15e3, 15e3]))
    else:
        b0 = window_bounds(area, UTM)
    base, bb = read(P / "dem_gem15_dtm_utm.vrt", b0)
    hs = hillshade(base, (bb[2] - bb[0]) / base.shape[1])
    made = []
    for name, (f, label, cmap, rule, cap) in LAYERS.items():
        if only and name not in only:
            continue
        if name == "dsm_minus_dtm":
            a1, _ = read(P / "dem_gem15_dsm_utm.vrt", bb); a, _ = read(P / "dem_gem15_dtm_utm.vrt", bb); a = a1 - a
        elif name == "gem15_minus_glo30":
            a1, _ = read(P / "dem_gem15_dsm_utm.vrt", bb); g, _ = read(P / "dem_glo30_utm_full.tif", bb, maxpx=max(a1.shape))
            from scipy.ndimage import zoom
            g = zoom(g, (a1.shape[0] / g.shape[0], a1.shape[1] / g.shape[1]), order=1); a = a1 - g
        else:
            if not f.exists():
                print("missing", f.name); continue
            a, _ = read(f, bb, maxpx=max(base.shape))
            if a.shape != base.shape:
                from scipy.ndimage import zoom
                a = zoom(np.nan_to_num(a, nan=np.nanmedian(a)), (base.shape[0] / a.shape[0], base.shape[1] / a.shape[1]), order=1)
        if rule == "log":
            a = np.log10(np.where(a > 0, a, np.nan))
        lo, hi = rng(a, rule)
        asp = (bb[3] - bb[1]) / (bb[2] - bb[0])
        fig = plt.figure(figsize=(10, 10 * asp * 0.82 + 1.6))
        ax = fig.add_axes([0.03, 0.14, 0.82, 0.78])
        ext = [bb[0], bb[2], bb[1], bb[3]]
        ax.imshow(hs, cmap="gray", extent=ext, vmin=0, vmax=1, zorder=1)
        im = ax.imshow(a, cmap=cmap, vmin=lo, vmax=hi, extent=ext, alpha=0.62 if name != "elevation" else 0.55, zorder=2)
        overlays(ax, bb, vec); furniture(ax, bb)
        cax = fig.add_axes([0.87, 0.3, 0.02, 0.5]); cb = fig.colorbar(im, cax=cax); cb.set_label(label)
        title = AREAS[area] and f"{area.title()} ({AREAS[area][2] * 2} km across)" or "Karakoram Highway corridor"
        fig.text(0.03, 0.965, f"{label.split(' (')[0]}: {title}", fontsize=13, fontweight="bold")
        fig.text(0.03, 0.105, cap, fontsize=8, color="0.25", wrap=True)
        legend(fig)
        fig.text(0.99, 0.002, "UTM 43N. Hillshade from GEM-15 dtm (Google, CC-BY 4.0). Roads, dams: OpenStreetMap. Faults: GEM. Earthquakes: USGS.",
                 fontsize=6.5, color="0.5", ha="right")
        o = FIG / f"{name}_{area}.png"; fig.savefig(o, dpi=150); plt.close(fig); made.append(o); print(o.name, flush=True)
    # contact sheet
    if made:
        n = len(made); c = 3; r = int(np.ceil(n / c))
        fig, axs = plt.subplots(r, c, figsize=(15, 5.2 * r))
        for ax in axs.flat:
            ax.axis("off")
        for ax, m in zip(axs.flat, made):
            ax.imshow(plt.imread(m)); ax.set_title(m.stem, fontsize=9)
        fig.tight_layout(); fig.savefig(FIG / f"index_{area}.png", dpi=90); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--area", default="all"); ap.add_argument("--only", default="")
    a = ap.parse_args()
    for ar in (AREAS if a.area == "every" else [a.area]):
        main(ar, set(filter(None, a.only.split(","))))
