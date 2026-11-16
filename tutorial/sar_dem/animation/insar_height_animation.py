"""How radar phase becomes height: a reproducible, procedural animation, from geometry to a height read at any point.

    python tutorial/sar_dem/animation/insar_height_animation.py [--lon 74.82 --lat 36.32 --km 12 --fps 12 --seed 7]
                                                                [--points "74.80,36.31;74.86,36.34"]

Everything is computed, nothing is drawn by hand. Terrain is real: a Copernicus GLO-30 patch (default Attabad, Hunza,
from data/raw/dem/glo30, made by src/02_dem.py). The radar is SIMULATED on that terrain with the geometry of the real
Sentinel-1 pair used in the tutorial (C band, lambda 5.55 cm, B_perp 440 m, slant range ~850 km, look angle 39 deg),
so every stage can be compared with the truth. Same seed, same frames.

Scenes (each a run of frames):
  1  geometry: two passes, the baseline, two ranges to one point; phase = 4 pi (R1 - R2) / lambda
  2  pre-processing: two single-look complex images; each one's phase looks random (speckle), amplitude shows slopes
  3  interferogram: image 1 x conj(image 2); the random part cancels; flat-earth fringes dominate
  4  flat-earth removal with the orbits: the fringes left follow the terrain; one fringe = height of ambiguity
  5  coherence and filtering: steep slopes (layover, shadow) and noise; multilook; mask below 0.4
  6  unwrapping: quality-guided region growing from a reference pixel, adding whole cycles (2 pi k)
  7  phase to height: z = -lambda R0 sin(theta) / (4 pi B_perp) * phi; reconstructed DEM against the truth
  8  read any point: wrapped phase, cycles k, unwrapped phase, height, true height, error
Writes tutorial/sar_dem/animation/out/insar_height.mp4, .gif and a PNG per scene (for the manual), and a JSON record
of every parameter and the point readouts.
"""
import argparse
import heapq
import json
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import uniform_filter, zoom

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
R = HERE.parents[2]
OUT = HERE / "out"
LAM, BPERP, R0, THETA = 0.05546576, 440.0, 850e3, np.radians(39.0)
INK, PAPER, ACC = "#1d2a33", "#fbfaf7", "#c0392b"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "figure.facecolor": PAPER, "axes.facecolor": PAPER,
                     "axes.edgecolor": "0.6", "savefig.facecolor": PAPER})
CYC = plt.get_cmap("twilight_shifted")
TER = plt.get_cmap("terrain")


# ------------------------------------------------------------------ data and forward model
def load_dem(lon, lat, km, n=360):
    tile = R / f"data/raw/dem/glo30/Copernicus_DSM_COG_10_N{int(lat):02d}_00_E{int(lon):03d}_00_DEM.tif"
    dlat = km / 111.32 / 2; dlon = dlat / np.cos(np.radians(lat))
    with rasterio.open(tile) as r:
        a = r.read(1, window=from_bounds(lon - dlon, lat - dlat, lon + dlon, lat + dlat, r.transform)).astype("float64")
    z = zoom(a, (n / a.shape[0], n / a.shape[1]), order=1)
    return z, km * 1000 / n, (lon - dlon, lat - dlat, lon + dlon, lat + dlat)


def simulate(z, px, seed):
    rng = np.random.default_rng(seed)
    ny, nx = z.shape
    x = np.arange(nx) * px                                     # ground range (m), radar looks along +x
    gy, gx = np.gradient(z, px)
    local = THETA - np.arctan(gx)                              # local incidence in the range plane
    layover, shadow = local <= 0.02, local >= np.pi / 2 - 0.02
    amp = np.clip(np.cos(local), 0.02, 1) ** 2 * (1 + 0.15 * rng.standard_normal(z.shape))
    amp[layover] = 1.5; amp[shadow] = 0.02
    k = 4 * np.pi / LAM
    phi_flat = -k * BPERP * x[None, :] / (R0 * np.tan(THETA)) * np.ones((ny, 1))
    phi_topo = -k * BPERP * z / (R0 * np.sin(THETA))
    coh = np.clip(0.92 - 0.9 * (np.abs(np.arctan(gx)) / np.radians(55)) ** 2, 0.05, 0.95)
    coh[layover | shadow] = 0.05
    speckle = rng.uniform(-np.pi, np.pi, z.shape)               # common to both images: cancels in the interferogram
    noise = rng.standard_normal(z.shape) * np.sqrt(1 - coh ** 2) / coh * 0.9
    s1 = np.sqrt(amp) * np.exp(1j * speckle)
    s2 = np.sqrt(amp) * np.exp(1j * (speckle - phi_flat - phi_topo - noise))
    return dict(x=x, amp=amp, layover=layover, shadow=shadow, coh=coh, s1=s1, s2=s2, phi_flat=phi_flat, phi_topo=phi_topo,
                h_a=LAM * R0 * np.sin(THETA) / (2 * BPERP))


def fringe_coherence(ifg, s1, s2, L=5):
    """Coherence and multilooked phase in an L x L window AFTER removing the local fringe frequency (the phase slope
    estimated from neighbour products), as real processors do on steep terrain; plain boxcar averaging across dense
    fringes would make good pixels look incoherent."""
    sm = lambda a: uniform_filter(a.real, 7) + 1j * uniform_filter(a.imag, 7)
    gx = np.angle(sm(ifg[:, 1:] * np.conj(ifg[:, :-1]))); gx = np.pad(gx, ((0, 0), (0, 1)), mode="edge")
    gy = np.angle(sm(ifg[1:, :] * np.conj(ifg[:-1, :]))); gy = np.pad(gy, ((0, 1), (0, 0)), mode="edge")
    h = L // 2
    num = np.zeros(ifg.shape, complex); p1 = np.zeros(ifg.shape); p2 = np.zeros(ifg.shape)
    for dy in range(-h, h + 1):
        for dx in range(-h, h + 1):
            sh = lambda a: np.roll(np.roll(a, -dy, 0), -dx, 1)
            num += sh(ifg) * np.exp(-1j * (gx * dx + gy * dy))
            p1 += sh(np.abs(s1) ** 2); p2 += sh(np.abs(s2) ** 2)
    return num, np.abs(num) / np.sqrt(p1 * p2)


def wrap(p):
    return (p + np.pi) % (2 * np.pi) - np.pi


def unwrap_quality(w, q, mask, ref):
    """Quality-guided unwrapping: grow from ref, always next to the best-quality pixel, adding the whole cycles that keep
    each step under half a cycle. Returns unwrapped phase and the visiting order (for the animation)."""
    ny, nx = w.shape
    u = np.full(w.shape, np.nan); done = np.zeros(w.shape, bool)
    u[ref] = w[ref]; done[ref] = True
    heap, order = [], [ref]
    def push(i, j, pi, pj):
        for a, b in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            if 0 <= a < ny and 0 <= b < nx and mask[a, b] and not done[a, b]:
                heapq.heappush(heap, (-q[a, b], a, b, i, j))
    push(*ref, *ref)
    while heap:
        _, a, b, i, j = heapq.heappop(heap)
        if done[a, b]:
            continue
        u[a, b] = u[i, j] + wrap(w[a, b] - w[i, j]); done[a, b] = True; order.append((a, b))
        push(a, b, i, j)
    return u, order


# ------------------------------------------------------------------ drawing helpers
class Film:
    def __init__(self, d):
        self.d, self.n = d, 0
        shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)

    def shot(self, fig, hold=1):
        for _ in range(hold):
            fig.savefig(self.d / f"f{self.n:05d}.png", dpi=100); self.n += 1


def frame(title, step, caption):
    fig = plt.figure(figsize=(12.8, 7.2))
    fig.text(0.03, 0.95, f"{step}", fontsize=11, color=ACC, fontweight="bold", va="top")
    fig.text(0.03, 0.91, title, fontsize=18, color=INK, fontweight="bold", va="top")
    fig.text(0.03, 0.045, caption, fontsize=10.5, color="0.25", va="bottom", wrap=True)
    fig.text(0.97, 0.015, "Simulated Sentinel-1 pair (lambda 5.55 cm, B_perp 440 m) on real Copernicus GLO-30 terrain (resampled), Attabad, Hunza",
             fontsize=7.5, color="0.5", ha="right")
    return fig


def img(ax, a, cmap, vmin=None, vmax=None, title=""):
    ax.imshow(a, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=11, color=INK)


# ------------------------------------------------------------------ scenes
def scene_geometry(film, z, sim, fps):
    row = z[z.shape[0] // 2]; x = sim["x"] / 1000
    for t in np.linspace(0, 1, int(fps * 4)):
        fig = frame("Two passes see one point from slightly different places", "1  Geometry",
                    "The satellite passes twice, 12 days apart, on orbits B_perp = 440 m apart. The distance to a point differs between "
                    "passes by a fraction of a wavelength; the radar records that fraction as phase: phi = 4 pi (R1 - R2) / lambda. "
                    "Higher ground changes that difference, which is why phase carries height.")
        ax = fig.add_axes([0.06, 0.16, 0.88, 0.68])
        ax.fill_between(x, row.min() - 300, row, color="#c9b38f"); ax.plot(x, row, color=INK, lw=1)
        sx, sy = -4.0, row.max() + 2600
        b = 0.7
        p = int(len(x) * (0.25 + 0.5 * t)); px_, pz = x[p], row[p]
        for k, (dx, dy, c) in enumerate(((0, 0, INK), (b * 0.6, b * 0.8 * 400, ACC))):
            ax.plot(sx + dx, sy + dy, marker=(3, 0, -90), ms=16, color=c)
            ax.plot([sx + dx, px_], [sy + dy, pz], color=c, lw=1.2, alpha=0.8)
            ax.text(sx + dx - 0.25, sy + dy + 120, f"pass {k + 1}", color=c, ha="right")
        ax.annotate("", xy=(sx + b * 0.6, sy + b * 0.8 * 400), xytext=(sx, sy), arrowprops=dict(arrowstyle="<->", color="0.4"))
        ax.text(sx + 0.5, sy - 260, "baseline", color="0.4", fontsize=9)
        ax.plot(px_, pz, "o", color=ACC, ms=7)
        ax.text(px_, pz + 180, f"P: {pz:,.0f} m", ha="center", color=INK, fontsize=10, fontweight="bold")
        ax.set_xlim(-6, x[-1] + 0.5); ax.set_ylim(row.min() - 300, sy + 900)
        ax.set_xlabel("ground range (km)"); ax.set_ylabel("height (m)")
        ax.spines[["top", "right"]].set_visible(False)
        film.shot(fig); plt.close(fig)


def scene_slc(film, sim, fps):
    a1, a2 = np.abs(sim["s1"]) ** 2, np.abs(sim["s2"]) ** 2
    for k in range(int(fps * 4)):
        fig = frame("Pre-processing: two single-look complex images", "2  Pre-processing",
                    "After focusing, co-registration (aligning image 2 to image 1 to a fraction of a pixel) and burst stitching, each "
                    "pixel is a complex number. Amplitude (left) shows slopes: bright facing the radar, black in shadow. "
                    "Phase (right) of ONE image looks random: each pixel holds many scatterers. Alone, it is useless.")
        axs = [fig.add_axes([0.04 + i * 0.24, 0.17, 0.22, 0.65]) for i in range(4)]
        img(axs[0], 10 * np.log10(a1), "gray", -15, 3, "image 1 amplitude (dB)")
        img(axs[1], np.angle(sim["s1"]), CYC, -np.pi, np.pi, "image 1 phase")
        img(axs[2], 10 * np.log10(a2), "gray", -15, 3, "image 2 amplitude (dB)")
        img(axs[3], np.angle(sim["s2"]), CYC, -np.pi, np.pi, "image 2 phase")
        film.shot(fig); plt.close(fig)


def scene_ifg(film, sim, fps):
    ifg = sim["s1"] * np.conj(sim["s2"])
    ny = ifg.shape[0]
    for t in np.linspace(0, 1, int(fps * 4)):
        fig = frame("The interferogram: image 1 x conjugate of image 2", "3  Interferogram",
                    "Multiplying by the complex conjugate subtracts the phases. The random part, the same in both images, cancels. "
                    "What is left is dominated by the flat-earth ramp: dense fringes caused by the curvature of the reference "
                    "ellipsoid and the change of look angle across the swath.")
        ax = fig.add_axes([0.05, 0.17, 0.4, 0.68]); ax2 = fig.add_axes([0.55, 0.17, 0.4, 0.68])
        n = max(1, int(ny * t))
        show = np.full(ifg.shape, np.nan); show[:n] = np.angle(ifg[:n])
        img(ax, np.angle(sim["s1"]), CYC, -np.pi, np.pi, "one image's phase (random)")
        img(ax2, show, CYC, -np.pi, np.pi, "interferogram phase (wrapped, -pi..pi)")
        film.shot(fig); plt.close(fig)
    return ifg


def scene_flat(film, sim, ifg, fps):
    flat = np.exp(1j * sim["phi_flat"])
    out = ifg * flat
    for t in np.linspace(0, 1, int(fps * 4)):
        mix = ifg * np.exp(1j * sim["phi_flat"] * t)
        fig = frame("Remove the flat-earth phase, using the orbits", "4  Flat-earth removal",
                    f"The orbits are known to centimetres, so the ramp can be computed and removed. What remains follows the terrain: "
                    f"each fringe (one full colour cycle) is one height of ambiguity, h_a = lambda R0 sin(theta) / (2 B_perp) = "
                    f"{sim['h_a']:.0f} m. Fringes crowd on steep slopes; if a pixel spans more than half a fringe the phase aliases and cannot be unwrapped, which is why fine pixels (here ~14 m) matter in the Karakoram.")
        ax = fig.add_axes([0.05, 0.17, 0.4, 0.68]); ax2 = fig.add_axes([0.55, 0.17, 0.4, 0.68])
        img(ax, np.angle(mix), CYC, -np.pi, np.pi, f"ramp removed: {100 * t:.0f} %")
        img(ax2, np.angle(out) if t > 0.98 else np.angle(mix), CYC, -np.pi, np.pi, "topographic fringes")
        film.shot(fig); plt.close(fig)
    return out


def scene_coh(film, sim, ifg_t, fps):
    num, coh = fringe_coherence(ifg_t, sim["s1"], sim["s2"])
    filt = np.angle(num)
    mask = coh >= 0.4
    for t in np.linspace(0, 1, int(fps * 4)):
        fig = frame("Coherence, filtering and the mask", "5  Quality",
                    "Coherence (0 to 1) says how alike the two images are in a small window. It drops with noise, vegetation, snow, "
                    "and to nothing in layover and shadow. Averaging neighbours (5 x 5, after removing the local fringe slope so "
                    "steep terrain is not mistaken for noise) cleans the phase. Pixels below 0.4 are masked: their phase cannot be trusted, so it is not guessed.")
        axs = [fig.add_axes([0.04 + i * 0.32, 0.17, 0.29, 0.66]) for i in range(3)]
        img(axs[0], coh, "magma", 0, 1, "estimated coherence")
        img(axs[1], np.where(t > 0.5, filt, np.angle(ifg_t)), CYC, -np.pi, np.pi, "filtered phase" if t > 0.5 else "raw phase")
        img(axs[2], np.where(mask, filt, np.nan), CYC, -np.pi, np.pi, f"masked: {100 * mask.mean():.0f} % kept")
        film.shot(fig); plt.close(fig)
    return filt, coh, mask


def scene_unwrap(film, sim, ifg_t, filt, coh, mask, fps, z):
    """Two routes. Naive: unwrap the full topographic phase (fails on steep slopes: cycle slips). Standard: subtract the
    phase predicted by an existing, coarser DEM, unwrap only the small residual, then add the DEM back (what HyP3, ISCE
    and SNAP do)."""
    scale = -LAM * R0 * np.sin(THETA) / (4 * np.pi * BPERP)
    ref = np.unravel_index(np.argmax(np.where(mask, coh - 1e-4 * z, -1)), z.shape)
    u_naive, _ = unwrap_quality(filt, coh, mask, ref)
    err_naive = scale * (u_naive - u_naive[ref]) + z[ref] - z
    for _ in range(int(fps * 3)):
        fig = frame("Why not unwrap everything directly?", "6  Unwrapping",
                    "On 30-45 degree Karakoram slopes the terrain phase changes about 1 radian per pixel; with noise, a step "
                    "sometimes looks like less than half a cycle when it is more. One wrong whole cycle is a jump of one height "
                    f"of ambiguity ({sim['h_a']:.0f} m), and it propagates downstream: the right panel shows these cycle slips.")
        ax = fig.add_axes([0.05, 0.17, 0.4, 0.68]); ax2 = fig.add_axes([0.55, 0.17, 0.4, 0.68])
        img(ax, np.where(mask, filt, np.nan), CYC, -np.pi, np.pi, "wrapped terrain phase (dense)")
        img(ax2, err_naive, "RdBu_r", -300, 300, "naive result minus truth (m)")
        film.shot(fig); plt.close(fig)
    from scipy.ndimage import gaussian_filter
    zref = gaussian_filter(z, 6) + 4.0               # an older, coarser DEM (~80 m smoothing, 4 m bias)
    phi_ref = -4 * np.pi / LAM * BPERP * zref / (R0 * np.sin(THETA))
    res = ifg_t * np.exp(-1j * phi_ref)
    num, coh_r = fringe_coherence(res, sim["s1"], sim["s2"])
    fr = np.angle(num); m = coh_r >= 0.4
    u, order = unwrap_quality(fr, coh_r, m, ref)
    cuts = np.linspace(1, len(order), int(fps * 6)).astype(int)
    for c in cuts:
        part = np.full(z.shape, np.nan)
        idx = np.array(order[:c]); part[idx[:, 0], idx[:, 1]] = u[idx[:, 0], idx[:, 1]]
        fig = frame("Unwrap only what the old DEM did not know", "6  Unwrapping",
                    "Subtract the phase an existing coarse DEM predicts; the residual has few, wide fringes. Then grow from a "
                    "reliable reference pixel to the best-quality neighbour first, adding whole cycles (2 pi k) so each step "
                    "changes by less than half a cycle. SNAPHU solves the same problem globally (MCF start, TOPO cost mode).")
        axs = [fig.add_axes([0.04 + i * 0.32, 0.17, 0.29, 0.66]) for i in range(3)]
        img(axs[0], np.angle(np.exp(1j * phi_ref)), CYC, -np.pi, np.pi, "phase predicted by the old DEM")
        img(axs[1], np.where(m, fr, np.nan), CYC, -np.pi, np.pi, "residual (wrapped)")
        axs[1].plot(ref[1], ref[0], "*", color="white", ms=14, mec=INK)
        img(axs[2], part, "viridis", np.nanpercentile(u, 1), np.nanpercentile(u, 99), f"unwrapped: {c:,} of {len(order):,} pixels")
        film.shot(fig); plt.close(fig)
    return u, fr, ref, zref, scale


def scene_height(film, u, ref, z, zref, scale, fps):
    zh = zref + scale * (u - u[ref]) + (z[ref] - zref[ref])     # one control point fixes the level
    err = zh - z; e0 = zref - z
    nm = lambda e: 1.4826 * np.nanmedian(np.abs(e - np.nanmedian(e)))
    for t in np.linspace(0, 1, int(fps * 5)):
        fig = frame("From phase to height", "7  Phase to height",
                    f"Height = old DEM + {scale:.2f} m per radian x unwrapped residual; one control point fixes the level. Then "
                    f"geocoding maps radar geometry to the map. Error of the old DEM: NMAD {nm(e0):.1f} m; after adding the radar "
                    f"residual: median {np.nanmedian(err):+.1f} m, NMAD {nm(err):.1f} m, on {100 * np.isfinite(zh).mean():.0f} % of pixels.")
        axs = [fig.add_axes([0.04 + i * 0.32, 0.17, 0.29, 0.66]) for i in range(3)]
        lo, hi = np.nanpercentile(z, 1), np.nanpercentile(z, 99)
        img(axs[0], zref - z, "RdBu_r", -40, 40, "old DEM minus truth (m)")
        img(axs[1], np.where(np.isfinite(zh), zref + (zh - zref) * t, np.nan), TER, lo, hi, "radar height")
        img(axs[2], err, "RdBu_r", -40, 40, "radar height minus truth (m)")
        film.shot(fig); plt.close(fig)
    return zh


def scene_points(film, z, zh, filt, u, ref, scale, box, points, fps):
    ny, nx = z.shape
    w, s, e, n = box
    pix = [(int((n - la) / (n - s) * ny), int((lo - w) / (e - w) * nx)) for lo, la in points]
    pix = [(min(max(i, 0), ny - 1), min(max(j, 0), nx - 1)) for i, j in pix]
    path = []
    for a, b in zip(pix, pix[1:] + pix[:1]):
        for t in np.linspace(0, 1, int(fps * 2.5)):
            path.append((int(a[0] + (b[0] - a[0]) * t), int(a[1] + (b[1] - a[1]) * t)))
    rec = []
    for i, j in path:
        ok = np.isfinite(u[i, j])
        k = int(np.round((u[i, j] - filt[i, j]) / (2 * np.pi))) if ok else None
        fig = frame("Read the height at any point", "8  Any point",
                    "At each pixel: the wrapped residual phase, the whole cycles k found by unwrapping, the unwrapped phase, and "
                    "the height it gives. Masked pixels (layover, shadow, low coherence) give no height: the honest answer is "
                    "'no data', to be filled from the other pass direction or another DEM.")
        ax = fig.add_axes([0.04, 0.15, 0.5, 0.72])
        img(ax, zh, TER, np.nanpercentile(z, 1), np.nanpercentile(z, 99), "")
        ax.axhline(i, color=ACC, lw=0.8); ax.axvline(j, color=ACC, lw=0.8); ax.plot(j, i, "o", mfc="none", mec=ACC, ms=14, mew=2)
        lon = w + (j + 0.5) / nx * (e - w); lat = n - (i + 0.5) / ny * (n - s)
        lines = [("position", f"{lat:.4f} N  {lon:.4f} E"), ("residual phase", f"{filt[i, j]:+.2f} rad")]
        if ok:
            lines += [("whole cycles k", f"{k:+d}"), ("unwrapped phase", f"{u[i, j]:+.2f} rad"),
                      ("relative to ref.", f"{u[i, j] - u[ref]:+.2f} rad"), ("x m per rad", f"{scale:.2f}"),
                      ("radar height", f"{zh[i, j]:,.0f} m"), ("true height", f"{z[i, j]:,.0f} m"), ("error", f"{zh[i, j] - z[i, j]:+.1f} m")]
        else:
            lines += [("whole cycles k", "-"), ("radar height", "no data (masked)"), ("true height", f"{z[i, j]:,.0f} m")]
        for r_, (a, b) in enumerate(lines):
            fig.text(0.59, 0.8 - r_ * 0.065, a, fontsize=12, color="0.4")
            fig.text(0.76, 0.8 - r_ * 0.065, b, fontsize=13, color=INK, fontweight="bold")
        film.shot(fig); plt.close(fig)
        if (i, j) in pix:
            rec.append(dict(lines))
    return rec


def main(a):
    OUT.mkdir(parents=True, exist_ok=True)
    z, px, box = load_dem(a.lon, a.lat, a.km)
    sim = simulate(z, px, a.seed)
    film = Film(OUT / "frames")
    marks = {}
    def mark(name):
        marks[name] = film.n
    mark("1_geometry"); scene_geometry(film, z, sim, a.fps)
    mark("2_slc"); scene_slc(film, sim, a.fps)
    mark("3_interferogram"); ifg = scene_ifg(film, sim, a.fps)
    mark("4_flat_earth"); ifg_t = scene_flat(film, sim, ifg, a.fps)
    mark("5_coherence"); filt, coh, mask = scene_coh(film, sim, ifg_t, a.fps)
    mark("6_unwrapping"); u, filt, ref, zref, scale = scene_unwrap(film, sim, ifg_t, filt, coh, mask, a.fps, z)
    mark("7_height"); zh = scene_height(film, u, ref, z, zref, scale, a.fps)
    pts = [tuple(map(float, p.split(","))) for p in a.points.split(";")] if a.points else \
          [(box[0] + (box[2] - box[0]) * fx, box[1] + (box[3] - box[1]) * fy) for fx, fy in ((0.3, 0.3), (0.7, 0.4), (0.5, 0.75), (0.2, 0.7))]
    mark("8_points"); rec = scene_points(film, z, zh, filt, u, ref, scale, box, pts, a.fps)
    names = list(marks); ends = [marks[k] for k in names[1:]] + [film.n]
    for name, end in zip(names, ends):
        last = end - 1 if name != "8_points" else marks[name]
        shutil.copy2(OUT / "frames" / f"f{last:05d}.png", OUT / f"scene_{name}.png")
    ff = shutil.which("ffmpeg")
    subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", str(OUT / "frames/f%05d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", str(OUT / "insar_height.mp4")], check=True)
    subprocess.run([ff, "-y", "-loglevel", "error", "-i", str(OUT / "insar_height.mp4"), "-vf",
                    "fps=8,scale=800:-1:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse", str(OUT / "insar_height.gif")], check=True)
    (OUT / "record.json").write_text(json.dumps({"lon": a.lon, "lat": a.lat, "km": a.km, "pixel_m": px, "seed": a.seed,
        "lambda_m": LAM, "bperp_m": BPERP, "R0_m": R0, "look_deg": 39.0, "h_ambiguity_m": sim["h_a"], "m_per_rad": scale,
        "frames": film.n, "scenes": marks, "points": rec}, indent=1))
    shutil.rmtree(OUT / "frames")
    print("frames", film.n, "h_a", round(sim["h_a"], 1), "m/rad", round(scale, 2), "->", OUT / "insar_height.mp4")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lon", type=float, default=74.82); ap.add_argument("--lat", type=float, default=36.32)
    ap.add_argument("--km", type=float, default=5); ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--seed", type=int, default=7); ap.add_argument("--points", default="")
    main(ap.parse_args())
