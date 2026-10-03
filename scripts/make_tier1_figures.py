"""Generate all Tier 1 (verification benchmark) figures from saved results.

    python scripts/make_tier1_figures.py

Reads ``results/{baseline,convergence,domain}`` and writes ``figures/tier1_*.png``
(+ ``.txt`` captions).  Nothing here runs a simulation.
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from em_rf_lab import analytical as A
from em_rf_lab import convergence as C
from em_rf_lab.io import load_run
from em_rf_lab.plotting import apply_style, save_fig
from em_rf_lab.validation import compare_to_mie

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
RES = REPO / "results"


def _ladder(pol="TM", short=False):
    root = RES / "shortrun" / "ladder" if short else RES / "convergence"
    runs = [load_run(p) for p in glob.glob(str(root / f"ladder_{pol}_*"))]
    runs.sort(key=lambda r: r.config.resolution)
    return runs


def fig_geometry():
    from em_rf_lab.config import SimulationConfig

    c = SimulationConfig()
    H, ax_h = 7.8, 6.3
    ax_w = ax_h * (c.cell_x + 2) / (c.cell_y + 2)            # equal aspect: narrow axes, centred on the common centre line
    fig = plt.figure(figsize=(CANVAS_W, H))
    ax = fig.add_axes([(CX - ax_w / 2) / CANVAS_W, 0.95 / H, ax_w / CANVAS_W, ax_h / H])
    ax.add_patch(plt.Rectangle((-c.cell_x / 2, -c.cell_y / 2), c.cell_x, c.cell_y, fc="#f4f4f4", ec="k"))
    inner = plt.Rectangle((-c.cell_x / 2 + c.pml, -c.cell_y / 2 + c.pml), c.cell_x - 2 * c.pml, c.cell_y - 2 * c.pml, fc="white", ec="#888", ls="--")
    ax.add_patch(inner)
    ax.add_patch(plt.Circle((0, 0), c.radius, fc="#9ecae1", ec="#08519c", label=f"cylinder (R={c.radius:g}, eps_r={c.eps_r:g})"))
    h = c.n2f_half_width
    ax.add_patch(plt.Rectangle((-h, -h), 2 * h, 2 * h, fill=False, ec="#d62728", lw=1.5, label="near-to-far / flux box"))
    ax.plot([c.source_x] * 2, [-c.cell_y / 2, c.cell_y / 2], color="#2ca02c", lw=2, label=f"line source (Gaussian profile, w={c.beam_waist:g})")
    ax.annotate("", xy=(c.source_x + 3.5, 0), xytext=(c.source_x + 0.3, 0), arrowprops=dict(arrowstyle="->", color="#2ca02c"))
    ax.text(-c.cell_x / 2 + 0.3, c.cell_y / 2 - 2.3, "PML", color="#555")
    ax.set_xlim(-c.cell_x / 2 - 1, c.cell_x / 2 + 1)
    ax.set_ylim(-c.cell_y / 2 - 1, c.cell_y / 2 + 1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (units of cylinder radius)")
    ax.set_ylabel("y (units of cylinder radius)")
    ax.set_title("Computational domain (production configuration)")
    h_, l_ = ax.get_legend_handles_labels()
    ax.legend(h_, [_wrap(t, 26) for t in l_], loc="upper left", bbox_to_anchor=(1.06, 0.97), fontsize=14, borderaxespad=0.0)
    save_fig(fig, FIG / "tier1_geometry.png", "Computational domain. The source line, closed near-to-far/flux box and PML are shown to scale; the dashed rectangle is the PML-free region.")


def fig_field_map():
    paths = glob.glob(str(RES / "baseline" / "baseline_*"))
    if not paths:
        return
    r = load_run(paths[0])
    if "map_F" not in r.arrays:
        return
    x, y, F = r["map_x"], r["map_y"], r["map_F"]
    ka = float(r["map_ka"][0])
    X, Y = np.meshgrid(x, y, indexing="ij")
    ref = A.near_field(X, Y, ka, r.config.radius, r.config.eps_r, "TM")
    err = np.abs(F - ref)
    rel = np.sqrt(np.sum(np.abs(F - ref) ** 2) / np.sum(np.abs(ref) ** 2))
    fig, axs, H = stack([5.0] * 3, gap=1.1, hscale=1.0)
    ext = [x[0], x[-1], y[0], y[-1]]
    vmax = max(np.abs(F.real).max(), np.abs(ref.real).max())
    for a, d, t in ((axs[0], F.real, "FDTD  Re E_z"), (axs[1], ref.real, "Exact (Mie)  Re E_z")):
        im = a.imshow(d.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        a.add_patch(plt.Circle((0, 0), r.config.radius, fill=False, ec="k", lw=1))
        a.set_title(t); a.set_xlabel("x"); a.set_ylabel("y"); a.grid(False)
    y_p1 = axs[0].get_position().y1 * H; y_p2 = axs[1].get_position().y0 * H
    span = (y_p1 - y_p2) * 0.85; mid = (y_p1 + y_p2) / 2
    colorbar_right(fig, H, im, mid - span / 2, mid + span / 2, "E_z / E_inc(0)")
    im2 = axs[2].imshow(err.T, origin="lower", extent=ext, cmap="magma")
    axs[2].add_patch(plt.Circle((0, 0), r.config.radius, fill=False, ec="w", lw=1))
    axs[2].set_title(f"|FDTD - exact|   (rel. L2 = {rel:.2e})"); axs[2].set_xlabel("x"); axs[2].grid(False)
    y3b = axs[2].get_position().y0 * H; y3t = axs[2].get_position().y1 * H
    colorbar_right(fig, H, im2, y3b + 0.075 * (y3t - y3b), y3t - 0.075 * (y3t - y3b), "absolute difference")
    save_fig(fig, FIG / "tier1_field_map.png", f"Complex near field at ka={ka:g} (TM, eps_r={r.config.eps_r:g}, res {r.config.resolution:g}); FDTD DFT field normalised by the incident amplitude vs the exact interior/exterior series. Relative L2 difference over the patch: {rel:.3e}.")


def fig_patterns(run=None):
    run = run or _ladder()[-1]
    c = compare_to_mie(run)
    kas = [2.0, 3.0, 4.0, 5.5]
    idx = [int(np.argmin(abs(c.ka - k))) for k in kas]
    fig, axs, H = stack([3.4] * 4, gap=1.25, top=1.3)
    ang = np.degrees(c.angles)
    for a, i in zip(axs, idx):
        a.semilogy(ang, c.w_ref[i], "k-", lw=2, label="exact (Mie)")
        a.semilogy(ang, run["echo_width"][i], "r--", lw=1.6, label=f"FDTD (res {run.config.resolution:g})")
        a.set_title(f"ka = {c.ka[i]:.1f}   rel. L2 = {c.rel_l2_angle[i]:.1e}")
        a.set_xlabel("scattering angle (deg; 0 = forward)")
        a.set_xlim(0, 180)
    for a in axs:
        a.set_ylabel("echo width W (units of R)")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(CX / CANVAS_W, 1 - 0.15 / H), ncol=2, fontsize=14, frameon=True)
    wrap_labels(axs)
    save_fig(fig, FIG / "tier1_patterns.png", "Far-field echo width W(phi)=2*pi*rho*|F_s|^2/|F_i|^2 from the FDTD near-to-far transform against the exact Mie series, TM polarization, eps_r=4.")


def fig_convergence():
    from em_rf_lab.validation import fit_frequency_shift

    runs = _ladder()
    short = _ladder(short=True)
    res = np.array([r.config.resolution for r in runs])
    comps = [compare_to_mie(r) for r in runs]
    rms = np.array([c.rms_rel_l2_angle for c in comps])
    spec = np.array([c.rel_l2_spectrum for c in comps])
    sfx = np.array([np.max(np.abs(c.sigma_err_flux)) for c in comps])
    shifts = [fit_frequency_shift(r) for r in runs]
    delta = np.array([abs(s.delta) for s in shifts])
    rms_comp = np.array([s.rms_l2_compensated for s in shifts])
    res_s = np.array([r.config.resolution for r in short])
    rms_s = np.array([compare_to_mie(r).rms_rel_l2_angle for r in short])
    fit = C.fit_order(res[res >= 15], rms[res >= 15]) if (res >= 15).sum() >= 2 else None
    fit_d = C.fit_order(res, delta)
    fig, axs, H = stack([3.9] * 3, gap=1.3)
    a = axs[0]
    a.loglog(res, rms, "o-", label="time-converged runs (this repo's production setting)")
    a.loglog(res_s, rms_s, "s--", color="gray", label="short runs (decay 1e-5, interior probe): duration artifact")
    a.loglog(res, rms_comp, "^-", label="converged, after removing fitted rigid frequency shift")
    a.loglog(res, spec, "v:", label="converged: rel. L2 error of sigma_s(ka)")
    rr = np.array([res.min(), res.max()])
    a.loglog(rr, rms[0] * (res[0] / rr) ** 2, color="k", ls=":", lw=0.8, label="ideal h^2")
    a.set_xlabel("FDTD resolution (pixels per unit length)")
    a.set_ylabel("error vs exact Mie solution")
    a.set_title("Grid convergence (TM, eps_r=4, ka in [2,6])")
    legend_right(a, 24)
    a = axs[1]
    a.loglog(res, delta, "o-", label="|fitted shift| delta")
    a.loglog(rr, fit_d.predict(1 / rr), "k:", label=f"fit: delta ~ h^{fit_d.p:.2f}")
    a.set_xlabel("FDTD resolution")
    a.set_ylabel("fractional frequency shift of resonances")
    a.set_title("Numerical dispersion: rigid spectral shift")
    legend_right(a, 24)
    a = axs[2]
    diffs = []
    for r0, r1 in zip(runs[:-1], runs[1:]):
        num = np.sqrt(np.sum((r0["echo_width"] - r1["echo_width"]) ** 2)); den = np.sqrt(np.sum(r1["echo_width"] ** 2))
        diffs.append(num / den)
    a.loglog(res[:-1], diffs, "o-", label="eps_h = ||W_h - W_next|| / ||W_next||")
    a.set_xlabel("resolution (coarser member of the pair)")
    a.set_ylabel("successive-grid difference")
    a.set_title("Self-convergence (no analytical reference)")
    legend_right(a, 24)
    wrap_labels(axs)
    local = C.local_orders(res, rms)
    txt = "local orders of the converged RMS error: " + ", ".join(f"{p:.2f}" for p in local)
    save_fig(fig, FIG / "tier1_convergence.png", f"Grid convergence of FDTD vs the exact series. Fitted dispersion-shift order {fit_d.p:.2f}; {txt}.")
    return dict(res=res.tolist(), rms=rms.tolist(), spec=spec.tolist(), flux=sfx.tolist(), shift=[s.delta for s in shifts], rms_compensated=rms_comp.tolist(),
                shift_order=fit_d.p, local_orders=local.tolist(), self_diff=[float(d) for d in diffs],
                shortrun=dict(res=res_s.tolist(), rms=rms_s.tolist()),
                order_fit=(dict(p=fit.p, r2=fit.r2) if fit else None))


def fig_error_vs_runtime():
    runs = _ladder()
    wall = np.array([r.meta["wall_time_total_s"] for r in runs])
    res = np.array([r.config.resolution for r in runs])
    rms = np.array([compare_to_mie(r).rms_rel_l2_angle for r in runs])
    fig, (ax,), H = stack([4.2], top=0.8, bottom=1.4)
    ax.loglog(wall, rms, "o-")
    for w, e, rs in zip(wall, rms, res):
        ax.annotate(f"{rs:g}", (w, e), textcoords="offset points", xytext=(6, 5), fontsize=12)
    ax.set_xlabel("wall-clock time per configuration, serial, 1 core (s)")
    ax.set_ylabel("RMS rel. L2 error of W(phi)")
    ax.set_title("Accuracy-cost tradeoff (labels: resolution)")
    wrap_labels(ax)
    save_fig(fig, FIG / "tier1_error_vs_runtime.png", "Error versus wall-clock cost of one full reference+scattering configuration (41 frequencies, 91 angles, serial).")


def fig_spectrum():
    run = _ladder()[-1]
    c = compare_to_mie(run)
    fine = np.linspace(run.config.ka_min, run.config.ka_max, 801)
    fig, axs, H = stack([4.2, 2.0], gap=0.35, sharex=True)
    axs[0].tick_params(labelbottom=False)
    axs[0].plot(fine, [A.total_scattering_width(k, 1, run.config.eps_r, "TM") for k in fine], "k-", label="exact (Mie series)")
    axs[0].plot(run["ka"], run["sigma_s_n2f"], "ro", ms=4, label="FDTD (near-to-far)")
    axs[0].plot(run["ka"], run["sigma_s_flux"], "b+", ms=6, label="FDTD (flux)")
    axs[0].set_ylabel("total scattering width sigma_s (units of R)")
    axs[0].set_title(f"Scattering spectrum, TM, eps_r = {run.config.eps_r:g} (res {run.config.resolution:g})")
    legend_right(axs[0], 22)
    axs[1].plot(c.ka, 100 * c.sigma_err_n2f, "ro-", ms=3, label="near-to-far")
    axs[1].plot(c.ka, 100 * c.sigma_err_flux, "b+-", ms=4, label="flux")
    axs[1].axhline(0, color="k", lw=0.5)
    axs[1].set_xlabel("size parameter ka = 2 pi f R / c")
    axs[1].set_ylabel("error (%)")
    legend_right(axs[1], 22)
    wrap_labels(axs)
    save_fig(fig, FIG / "tier1_spectrum.png", "Total scattering width spectrum from two independent FDTD estimators (near-to-far integral and scattered-flux) against the exact series; lower panel: signed relative error.")


def fig_domain():
    """Sensitivity of the far-field result to each numerical-setup parameter.

    Metric: relative L2 difference of W(phi) (RMS over ka, and for ka <= 3) with respect to the
    production setting, i.e. how much the answer *moves*; meaningful even when grid error dominates.
    """
    paths = sorted(glob.glob(str(RES / "domain" / "dom_*")))
    if not paths:
        return {}
    runs = {load_run(p).config.name: load_run(p) for p in paths}
    base = runs["dom_pml5"]
    groups = {
        "pml": ("dom_pml", "PML thickness"),
        "beam_waist": ("dom_waist", "Gaussian source waist w"),
        "n2f_half_width": ("dom_n2f", "near-to-far box half-width"),
        "min_run_time": ("dom_T", "run time after source"),
        "source_fwidth_factor": ("dom_bw", "source bandwidth factor"),
    }

    def diff(r):
        ka = base["ka"]
        d = np.array([np.linalg.norm(r["echo_width"][i] - base["echo_width"][i]) / np.linalg.norm(base["echo_width"][i]) for i in range(len(ka))])
        return float(np.sqrt(np.mean(d**2))), float(np.sqrt(np.mean(d[ka <= 3.0] ** 2))), float(np.max(d))

    out = {}
    fig, axs, H = stack([3.2] * 5, gap=1.25, top=2.2, sharey=True)
    for a, (attr, (prefix, label)) in zip(axs, groups.items()):
        sel = sorted([r for n, r in runs.items() if n.startswith(prefix)], key=lambda r: getattr(r.config, attr))
        vals = [getattr(r.config, attr) for r in sel]
        d = [diff(r) for r in sel]
        a.semilogy(range(len(vals)), [max(x[0], 1e-5) for x in d], "s--", label="RMS over all ka")
        a.semilogy(range(len(vals)), [max(x[1], 1e-5) for x in d], "o-", label="RMS over ka <= 3")
        a.set_xticks(range(len(vals)))
        a.set_xticklabels([f"{v:g}" for v in vals])
        a.set_xlabel(label)
        out[attr] = dict(values=vals, rms_all=[x[0] for x in d], rms_lowka=[x[1] for x in d], max=[x[2] for x in d])
    for a in axs[:len(groups)]:
        a.set_ylabel("change in W(phi) vs production setting (rel. L2)")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.suptitle(_wrap("Numerical-setup sensitivity at fixed resolution (points equal to the production value show ~0 by construction)", 62),
                 x=CX / CANVAS_W, y=1 - 0.1 / H, va="top", fontsize=14)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(CX / CANVAS_W, 1 - 1.0 / H), ncol=2, fontsize=14, frameon=True)
    wrap_labels(axs)
    save_fig(fig, FIG / "tier1_domain_studies.png", "One-factor-at-a-time sensitivity at fixed resolution: how much the far-field result changes when each numerical parameter is varied around the production setting.")
    return out


if __name__ == "__main__":
    apply_style()
    set_big_fonts()
    fig_geometry()
    fig_field_map()
    fig_patterns()
    summary = fig_convergence()
    fig_error_vs_runtime()
    fig_spectrum()
    summary["domain"] = fig_domain()
    (RES / "tier1_summary.json").write_text(json.dumps(summary, indent=2))
    print("figures written to", FIG)
