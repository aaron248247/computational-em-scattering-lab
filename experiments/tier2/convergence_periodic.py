"""Tier 2 verification: convergence of resonance metrics (f0, Q, peak R) for the cylinder array.

Two representative resonances (one per polarisation) are re-simulated while varying, one at a time,
the grid resolution, the structure-to-PML margin, the PML thickness and the R/T monitor gap.
Results: ``results/tier2/convergence`` (raw runs) and ``results/tier2/convergence_summary.json``.

    python experiments/tier2/convergence_periodic.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from em_rf_lab import resonance as Rz
from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "tier2" / "convergence"
CASES = {  # label: (polarisation, pitch, f window)
    "TM_L3": ("TM", 3.0, (0.18, 0.34)),
    "TE_L4": ("TE", 4.0, (0.14, 0.27)),
}
STUDIES = {
    "resolution": [20, 30, 40, 60],
    "margin": [4.0, 6.0, 9.0],
    "pml": [3.0, 4.0, 6.0],
    "monitor_gap": [2.0, 3.0, 5.0],
}
BASE = dict(resolution=40.0, margin=6.0, pml=4.0, monitor_gap=3.0)


def config(case: str, study: str, value: float) -> PeriodicConfig:
    pol, pitch, (f0, f1) = CASES[case]
    kw = dict(BASE)
    kw[study] = float(value)
    return PeriodicConfig(
        name=f"conv_{case}_{study}{value:g}", polarization=pol, pitch=pitch, f_min=f0, f_max=f1,
        n_freq=321, decay_tol=1e-7, max_run_time=2000.0, **kw,
    )


def metrics(run) -> dict:
    f, R, T = run["freq"], run["R"], run["T"]
    i = int(np.argmax(R))
    q_fwhm = Rz.fwhm_q(f, R, i)
    half = max(6 * (f[1] - f[0]), 1.5 * f[i] / q_fwhm) if np.isfinite(q_fwhm) else 0.02
    try:
        fit = Rz.fit_resonance(f, R, f_guess=float(f[i]), half_window=half, gamma_guess=half / 6)
        f_fit, q_fit, rms = fit.f0, fit.q, fit.rms_residual
    except Exception:
        f_fit = q_fit = rms = float("nan")
    return dict(f_peak=float(f[i]), r_peak=float(R[i]), q_fwhm=float(q_fwhm), f_fit=float(f_fit), q_fit=float(q_fit),
                fit_rms=float(rms), max_abs_A=float(np.max(np.abs(1 - R - T))), decay=bool(run.meta.get("decay_reached")))


if __name__ == "__main__":
    summary: dict = {}
    for case in CASES:
        for study, values in STUDIES.items():
            for v in values:
                r = run_periodic_or_load(config(case, study, v), str(OUT))
                m = metrics(r)
                summary.setdefault(case, {}).setdefault(study, {})[str(v)] = m
                print(f"{case} {study}={v:g}: f_fit={m['f_fit']:.5f} Q={m['q_fit']:.1f} Rpk={m['r_peak']:.4f} A={m['max_abs_A']:.1e}", flush=True)
    (REPO / "results" / "tier2" / "convergence_summary.json").write_text(json.dumps(summary, indent=1))
