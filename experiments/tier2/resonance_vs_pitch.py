"""Tier 2 primary study: lowest collective resonance versus pitch (fine frequency sampling).

The coarse survey in ``pitch_study.py`` locates the lowest strong reflection peak; here each is re-simulated
in a narrow window (321 points) so that f0, Q and the >=90 % reflection band are resolved by the line-shape fit.
Two independent Q estimates are produced: coupled-mode line-shape fit, and (``--harminv``) time-domain ringdown.
Results: ``results/tier2/lowest_resonance`` and ``results/tier2/lowest_resonance_summary.json``.

    python experiments/tier2/resonance_vs_pitch.py [--harminv]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from em_rf_lab import resonance as Rz
from em_rf_lab.io import load_run
from em_rf_lab.modes import find_modes
from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load
from em_rf_lab.tier2_analysis import band_edges

REPO = Path(__file__).resolve().parents[2]
COARSE = REPO / "results" / "tier2" / "pitch"
OUT = REPO / "results" / "tier2" / "lowest_resonance"
PITCHES = [2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]
# The coarse survey is unconverged for high-Q TM resonances at short pitch (it can select a spurious peak),
# so those window centres are taken by interpolating the converged neighbours (pitch 2.2 and 3.0).
CENTRE_OVERRIDE = {("TM", 2.4): 0.272, ("TM", 2.6): 0.267, ("TM", 2.8): 0.262}


def coarse_peak(pol: str, pitch: float) -> float:
    run = next(load_run(d, config_cls=PeriodicConfig) for d in sorted(COARSE.glob(f"pitch_{pol}_L{pitch:g}_*")))
    f, R = run["freq"], run["R"]
    idx = Rz.find_reflection_peaks(f, R, 0.3, 0.1)
    return float(f[int(idx[0])])


def fine_config(pol: str, pitch: float, f_c: float) -> PeriodicConfig:
    return PeriodicConfig(name=f"low_{pol}_L{pitch:g}", polarization=pol, pitch=pitch, f_min=0.75 * f_c, f_max=1.25 * f_c,
                          n_freq=321, decay_tol=1e-7, max_run_time=2500.0)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--harminv", action="store_true")
    ap.add_argument("--only", default="", help="comma list like TM_2.4,TM_2.6")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    summ = REPO / "results" / "tier2" / "lowest_resonance_summary.json"
    out: dict = json.loads(summ.read_text()) if summ.exists() else {}
    only = {x for x in args.only.split(",") if x}
    for pol in ("TM", "TE"):
        for p in PITCHES:
            if only and f"{pol}_{p:g}" not in only:
                continue
            fc = CENTRE_OVERRIDE.get((pol, p)) or coarse_peak(pol, p)
            cfg = fine_config(pol, p, fc)
            r = run_periodic_or_load(cfg, str(OUT))
            f, R, T = r["freq"], r["R"], r["T"]
            i = int(np.argmax(R))
            half = 0.2 * f[i]
            q_fwhm = Rz.fwhm_q(f, R, i)
            rec = dict(pol=pol, pitch=p, f_coarse=fc, f_peak=float(f[i]), r_peak=float(R[i]), q_fwhm=float(q_fwhm),
                       max_abs_A=float(np.max(np.abs(1 - R - T))), decay=bool(r.meta.get("decay_reached")))
            try:
                fit = Rz.fit_resonance(f, R, f_guess=float(f[i]), half_window=half, gamma_guess=half / 8)
                rec.update(f_fit=fit.f0, q_fit=fit.q, fit_rms=fit.rms_residual)
            except Exception as e:  # noqa: BLE001
                rec.update(f_fit=float("nan"), q_fit=float("nan"), fit_rms=float("nan"), fit_error=str(e))
            be = band_edges(f, R, i, 0.9)
            rec["band90_fractional"] = None if be is None else be.fractional
            if args.harminv:
                ms = find_modes(cfg, fcen=float(f[i]), df=0.3 * float(f[i]), run_time=1500.0, q_min=3.0)
                rec["harminv"] = [dict(freq=m.freq, q=m.q, amplitude=m.amplitude, error=m.error) for m in ms]
            out[f"{pol}_{p:g}"] = rec
            print(f"{pol} L={p:g}: f={rec['f_fit']:.4f} Q={rec['q_fit']:.1f} Rpk={rec['r_peak']:.4f} band90={rec['band90_fractional']}", flush=True)
            summ.write_text(json.dumps(out, indent=1))
