"""Generate the Tier 2 (periodic-array investigation) figures from saved results.

    python scripts/make_tier2_figures.py

Reads ``results/tier2/*`` and writes ``figures/tier2_*.png`` (+ ``.txt`` captions).  Runs no simulation.
"""

from __future__ import annotations

import glob
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from em_rf_lab import analytical as A
from em_rf_lab.io import load_run
from em_rf_lab.periodic import PeriodicConfig, slab_rt, wood_anomaly_frequencies
from em_rf_lab.plotting import apply_style, save_fig

# ---------------------------------------------------------------------------------------------
# Common figure layout.  Every figure has the same canvas width and the same horizontal axes box
# (left edge and width), so the plot areas line up when the PNGs are stacked in the documentation.
# Legends / colour bars live in the reserved right-hand margin; fonts are sized for a figure that is
# shrunk to document width.  Only layout is handled here: no data, labels or styles are changed.
# ---------------------------------------------------------------------------------------------
import textwrap
from matplotlib.ticker import LogFormatterSciNotation

CANVAS_W = 13.0                       # inches, identical for all figures
AX_LEFT = 1.7                         # left edge of the axes block (inches)
AX_W = 8.0                            # width of the axes block (inches)
HS = 1.2                              # vertical scale factor for stacked panel heights
CX = AX_LEFT + AX_W / 2               # shared horizontal centre line of the plots (inches)


def _wrap(text, n):
    return "\n".join(textwrap.wrap(str(text), n)) if text and len(str(text)) > n and "$" not in str(text) else text


def set_big_fonts():
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12, "legend.fontsize": 12, "figure.constrained_layout.use": False})


def stack(heights, gap=1.25, top=0.8, bottom=0.95, sharex=False, sharey=False, hscale=HS):
    """Vertical stack of axes (heights in inches) on the common horizontal axes box."""
    heights = [h * hscale for h in heights]
    n = len(heights)
    H = top + sum(heights) + (n - 1) * gap + bottom
    fig = plt.figure(figsize=(CANVAS_W, H))
    axs, y_top = [], H - top
    for i, h in enumerate(heights):
        y0 = y_top - h
        kw = {}
        if sharex and axs:
            kw["sharex"] = axs[0]
        if sharey and axs:
            kw["sharey"] = axs[0]
        axs.append(fig.add_axes([AX_LEFT / CANVAS_W, y0 / H, AX_W / CANVAS_W, h / H], **kw))
        y_top = y0 - gap
    return fig, axs, H


def legend_right(ax, wrap=24, **kw):
    """Legend in the right-hand margin, beside the axes (labels wrapped, not reworded)."""
    h, l = ax.get_legend_handles_labels()
    return ax.legend(h, [_wrap(t, wrap) for t in l], loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, **kw)


def wrap_labels(axs, ylab=30, xlab=44, title=46):
    for a in np.atleast_1d(axs).ravel():
        a.set_ylabel(_wrap(a.get_ylabel(), ylab))
        a.set_xlabel(_wrap(a.get_xlabel(), xlab))
        a.set_title(_wrap(a.get_title(), title))
        for axis, sc in ((a.xaxis, a.get_xscale()), (a.yaxis, a.get_yscale())):   # label minor ticks on short log axes
            if sc == "log":
                axis.set_minor_formatter(LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(2, 0.4)))


def colorbar_right(fig, H, im, y_bottom, y_top, label, width=0.3, pad=0.3):
    """Colour bar in the right margin spanning [y_bottom, y_top] (inches from the bottom of the canvas)."""
    cax = fig.add_axes([(AX_LEFT + AX_W + pad) / CANVAS_W, y_bottom / H, width / CANVAS_W, (y_top - y_bottom) / H])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(label)
    return cb


REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "figures"
T2 = REPO / "results" / "tier2"
PITCHES = [2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]
COL = {"TM": "#1f77b4", "TE": "#d62728"}


def _load(pattern):
    return [load_run(p, config_cls=PeriodicConfig) for p in sorted(glob.glob(str(T2 / pattern)))]


def fig_slab():
    runs = {r.config.name: r for r in _load("slab/*")}
    fig, ax, H = stack([4.0, 4.0, 4.0], gap=1.3)
    r = runs["slab_TM_normal_res40"]
    c = r.config
    f = r["freq"]
    Rr, Tr = slab_rt(f, math.sqrt(c.eps_r), c.thickness, 0.0, "TM")
    ax[0].plot(f, Rr, "k", lw=2, label="exact (Airy) R")
    ax[0].plot(f, r["R"], "r--", label="FDTD R (res 40)")
    ax[0].plot(f, Tr, "gray", lw=2, label="exact T")
    ax[0].plot(f, r["T"], "b--", label="FDTD T")
    ax[0].set(xlabel="frequency f (c/length)", ylabel="R, T", title="Homogeneous slab, normal incidence")
    legend_right(ax[0], 22)
    res, err = [], []
    for res_ in (20, 30, 40, 60, 80):
        r = runs[f"slab_TM_normal_res{res_}"]
        Rr, _ = slab_rt(r["freq"], math.sqrt(c.eps_r), c.thickness, 0.0, "TM")
        res.append(res_)
        err.append(float(np.sqrt(np.mean((r["R"] - Rr) ** 2))))
    p = np.polyfit(np.log(res), np.log(err), 1)
    ax[1].loglog(res, err, "o-", label="RMS error of R")
    ax[1].loglog(res, np.exp(p[1]) * np.array(res, float) ** p[0], "k:", label=f"fit: order {-p[0]:.2f}")
    ax[1].set(xlabel="resolution", ylabel="RMS |R_FDTD - R_exact|", title="Convergence of the periodic-cell workflow")
    legend_right(ax[1], 22)
    for pol, ls in (("TM", "-"), ("TE", "--")):
        for ky, col in ((0.5, "C0"), (1.0, "C1"), (1.5, "C2")):
            r = runs[f"slab_{pol}_ky{ky:g}_res40"]
            Rr, _ = slab_rt(r["freq"], math.sqrt(c.eps_r), c.thickness, ky, pol)
            ax[2].plot(r["freq"], r["R"] - Rr, ls, color=col, label=f"{pol} ky={ky:g}")
    ax[2].set(xlabel="frequency", ylabel="R_FDTD - R_exact", title="Oblique incidence (fixed ky)")
    legend_right(ax[2], 22, ncol=1)
    wrap_labels(ax)
    save_fig(fig, FIG / "tier2_slab_validation.png",
             "Verification of the periodic-cell reflectance/transmittance workflow against the exact slab (Airy) solution: "
             f"spectra (left), second-order grid convergence (centre, fitted order {-p[0]:.2f}) and oblique-incidence residuals (right).")


def fig_pitch_maps():
    fig, ax, H = stack([4.0, 4.0], gap=1.3, sharex=True)
    for a, pol in zip(ax, ("TM", "TE")):
        rows, freq = [], None
        for p in PITCHES:
            r = _load(f"pitch/pitch_{pol}_L{p:g}_*")[0]
            freq = r["freq"]
            rows.append(r["R"])
        im = a.pcolormesh(freq, np.arange(len(PITCHES)), np.array(rows), vmin=0, vmax=1, cmap="viridis", shading="nearest")
        a.set_yticks(range(len(PITCHES)), [f"{p:g}" for p in PITCHES])
        for j, p in enumerate(PITCHES):
            a.plot([1 / p] * 2, [j - 0.5, j + 0.5], "w", lw=2)
        a.set(xlabel="frequency f (c/length)", title=f"{pol}: reflectance R(f), normal incidence")
    ax[0].set_ylabel("pitch Λ (cylinder radius = 1)")
    wrap_labels(ax)
    colorbar_right(fig, H, im, 0.95 + 0.2, H - 0.8 - 0.2, "R")
    save_fig(fig, FIG / "tier2_pitch_maps.png",
             "Coarse survey (551 frequencies) of normal-incidence reflectance of a cylinder array (R=1, eps_r=4) versus pitch. "
             "White ticks mark the first Rayleigh-Wood anomaly f=c/Λ. Under-sampled/unconverged high-Q lines at short TM pitch are "
             "re-examined in the finer studies.")


def fig_lowest_resonance():
    d = json.loads((T2 / "lowest_resonance_summary.json").read_text())
    hq = json.loads((T2 / "highq_summary.json").read_text())
    branch = {}
    for key, v in d.items():
        v = dict(v)
        pol, p = key.split("_")
        if pol == "TM" and p in hq and "f_ringdown" in hq[p]:      # short-pitch TM: narrow mode from long ringdown
            h = hq[p]
            v.update(f_fit=h["f_ringdown"], q_fit=h["q_ringdown"], q_fwhm=h.get("q_fwhm_narrow") or float("nan"),
                     decay=True, max_abs_A=0.0, r_peak=h["r_peak_narrow"], band90_fractional=None, harminv=[], source="ringdown")
        else:
            v["source"] = "line-shape fit"
        branch[key] = dict(pol=pol, pitch=float(p), f0=v["f_fit"], q=v["q_fit"], source=v["source"], r_peak=v["r_peak"],
                           band90_fractional=v["band90_fractional"])
        d[key] = v
    (T2 / "lowest_branch.json").write_text(json.dumps(branch, indent=1))
    fig, ax, H = stack([3.6, 3.6, 3.6], gap=1.3)
    for pol in ("TM", "TE"):
        pts = [(p, d[f"{pol}_{p:g}"]) for p in PITCHES if f"{pol}_{p:g}" in d]
        good = [(p, v) for p, v in pts if v["decay"] and v["max_abs_A"] < 0.01 and v["r_peak"] > 0.95]
        bad = [(p, v) for p, v in pts if (p, v) not in good]
        for sel, mk, fc in ((good, "o", COL[pol]), (bad, "o", "none")):
            if not sel:
                continue
            P = [p for p, _ in sel]
            ax[0].plot(P, [v["f_fit"] for _, v in sel], mk, color=COL[pol], mfc=fc, label=pol if fc != "none" else None)
            ax[1].semilogy(P, [v["q_fit"] for _, v in sel], mk, color=COL[pol], mfc=fc)
        ax[1].semilogy([p for p, v in pts if v["decay"] and np.isfinite(v["q_fwhm"])], [v["q_fwhm"] for p, v in pts if v["decay"] and np.isfinite(v["q_fwhm"])], "^", color=COL[pol], alpha=.45, ms=5)
        gp = [p for p, _ in good]
        ax[0].plot(gp, [v["f_fit"] for _, v in good], "-", color=COL[pol], lw=1, alpha=.6)
        ax[1].semilogy(gp, [v["q_fit"] for _, v in good], "-", color=COL[pol], lw=1, alpha=.6)
        for p, v in pts:                       # Harminv cross-check where it matches the reflection peak
            for m in v.get("harminv", []):
                if abs(m["freq"] - v["f_fit"]) < 0.01 and m["q"] < 1e4:
                    ax[1].semilogy(p, m["q"], "x", color="k", ms=7)
        pb = [(p, v["band90_fractional"]) for p, v in pts if v["band90_fractional"]]
        ax[2].semilogy([p for p, _ in pb], [b for _, b in pb], "o-", color=COL[pol], label=pol)
    ax[1].semilogy([], [], "kx", label="Harminv ringdown")
    ax[1].semilogy([], [], "^", color="gray", label="FWHM estimate")
    ax[0].set(xlabel="pitch Λ", ylabel="f0 (c/length)", title="Resonance frequency")
    ax[1].set(xlabel="pitch Λ", ylabel="Q", title="Quality factor")
    ax[2].set(xlabel="pitch Λ", ylabel="fractional ≥90 % reflection bandwidth", title="Usable reflection band")
    ax[0].legend(); ax[1].legend(); ax[2].legend()
    wrap_labels(ax)
    save_fig(fig, FIG / "tier2_lowest_resonance.png",
             "Lowest reflection resonance of the cylinder array versus pitch: frequency, coupled-mode-fit Q (crosses: Harminv ringdown "
             "where it identifies the same mode; triangles: model-free FWHM estimate, which is lower for broad Fano-like lines) and fractional 90 % reflection bandwidth. Short-pitch TM points (2.2-2.8) are from long Harminv ringdown (Q at 2.4 and 2.6 approximate: ringdown only marginally converged); points failing a quality criterion "
             "(not decayed, energy residual >1 % or R_peak<0.95) and are not used for conclusions.")


def fig_convergence():
    s = json.loads((T2 / "convergence_summary.json").read_text())
    fig, ax, H = stack([3.9, 3.9], gap=1.3, bottom=2.6)
    for case, col in zip(s, ("C0", "C3")):
        for k, (study, vals) in enumerate(s[case].items()):
            xs = list(vals)
            f = np.array([vals[x]["f_fit"] for x in xs])
            q = np.array([vals[x]["q_fit"] for x in xs])
            ax[0].semilogy(range(len(xs)), np.abs(f / f[len(f) // 2 if study != "resolution" else -1] - 1) + 1e-9, "o-", label=f"{case}: {study}")
            ax[1].plot(range(len(xs)), q / q[-1 if study == "resolution" else len(q) // 2], "o-", label=f"{case}: {study}")
    ax[0].set(xlabel="study point index (see summary json)", ylabel="|f0/f0,ref − 1|", title="Resonance-frequency sensitivity")
    ax[1].set(xlabel="study point index", ylabel="Q / Q_ref", title="Q sensitivity")
    handles, labels = ax[0].get_legend_handles_labels()
    wrap_labels(ax)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(CX / CANVAS_W, 0.15 / H), ncol=2, fontsize=12, frameon=True)
    save_fig(fig, FIG / "tier2_convergence.png",
             "Sensitivity of f0 and Q to grid resolution (20-60), structure-PML margin, PML thickness and monitor gap, for a TM (pitch 3) and a TE (pitch 4) resonance.")


def fig_bands(pitch=3.0, pol="TM"):
    runs = sorted(_load(f"bands/band_{pol}_L{pitch:g}_ky*"), key=lambda r: r.config.ky)
    if len(runs) < 3:
        return
    modes = json.loads((T2 / "bands" / f"modes_L{pitch:g}_res40.json").read_text())
    ky = np.array([r.config.ky for r in runs])
    freq = runs[0]["freq"]
    Rm = np.array([np.where(freq > r.config.ky / (2 * np.pi), r["R"], np.nan) for r in runs])   # no propagating plane wave below the light line
    fig, ax, H = stack([4.2, 4.2], gap=1.3, sharex=True)
    im = ax[0].pcolormesh(ky * pitch / np.pi, freq, Rm.T, vmin=0, vmax=1, cmap="viridis", shading="nearest")
    kk = np.linspace(0, math.pi / pitch, 200)
    ax[0].plot(kk * pitch / np.pi, kk / (2 * np.pi), "w--", lw=1, label="light line f = ky/2π")
    ax[0].set(xlabel="ky Λ / π", ylabel="frequency f (c/length)", title=f"{pol}, Λ={pitch:g}: reflectance R(f, ky) at fixed ky")
    ax[0].legend(loc="lower right")
    colorbar_right(fig, H, im, 0.95 + 4.2 * HS + 1.3, H - 0.8, "R")
    pts = [(r.config.ky * pitch / np.pi, m["freq"], m["q"]) for r in runs for m in modes.get(f"{pol}_{r.config.ky:.4f}", [])
           if 5 < m["q"] < 2000 and m["error"] < 1e-3]
    if pts:
        x, y, q = map(np.array, zip(*pts))
        sc = ax[1].scatter(x, y, c=np.log10(q), cmap="plasma", s=30)
        colorbar_right(fig, H, sc, 0.95, 0.95 + 4.2 * HS, "log10 Q (Harminv)")
    ax[1].plot(kk * pitch / np.pi, kk / (2 * np.pi), "k--", lw=1)
    ax[1].set(xlabel="ky Λ / π", title="Ringdown modes (Harminv; 5 < Q < 2000)")
    wrap_labels(ax)
    save_fig(fig, FIG / f"tier2_bands_{pol}_L{pitch:g}.png",
             "Exploratory: reflectance map versus Bloch wavenumber and Harminv ringdown modes for the TM array at pitch 3. "
             "Below the light line (f < ky/2π) no propagating plane wave exists, so R is masked. At fixed ky the incidence angle varies with frequency. Harminv returns some duplicated or spurious "
             "modes; only modes with fit error <1e-3 and 5<Q<2000 are plotted.")


if __name__ == "__main__":
    apply_style()
    set_big_fonts()
    for fn in (fig_slab, fig_pitch_maps, fig_lowest_resonance, fig_convergence):
        fn()
    fig_bands()
    print("figures written to", FIG)
