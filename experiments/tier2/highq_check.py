"""Tier 2: high-Q modes of the TM array at short pitch -- ringdown (Harminv) vs narrow-window reflectance.

At short pitch the lowest TM mode has Q of order 10^2-10^3, narrower than the 321-point survey grid, so its
reflection line is under-sampled there.  Each mode is (i) located by Harminv with a long ringdown and (ii) searched for in a
very finely sampled reflectance window.  A ringdown mode that is absent from the normal-incidence
reflectance is *dark* (decoupled from the zeroth-order plane wave by symmetry).
Results: ``results/tier2/highq`` and ``results/tier2/highq_summary.json``.

    python experiments/tier2/highq_check.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from em_rf_lab import resonance as Rz
from em_rf_lab.modes import find_modes
from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "tier2" / "highq"
GUESS = {2.2: 0.2777, 2.4: 0.2725, 2.6: 0.267, 2.8: 0.2611}   # from the Harminv scan in resonance_vs_pitch.py

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    res: dict = {}
    for pitch, fg in GUESS.items():
        base = PeriodicConfig(name=f"hq_TM_L{pitch:g}", polarization="TM", pitch=pitch)
        modes = find_modes(base, fcen=fg, df=0.02, run_time=4000.0, q_min=20.0, err_max=1e-4)
        rec = dict(harminv=[dict(freq=m.freq, q=m.q, amplitude=m.amplitude, error=m.error) for m in modes])
        if modes:
            m = max(modes, key=lambda m: m.amplitude)
            cfg = PeriodicConfig(name=f"hq_TM_L{pitch:g}_RT", polarization="TM", pitch=pitch, f_min=m.freq * 0.995, f_max=m.freq * 1.005,
                                 n_freq=161, decay_tol=1e-8, max_run_time=8000.0)
            r = run_periodic_or_load(cfg, str(OUT))
            f, R = r["freq"], r["R"]
            i = int(np.argmax(R))
            rec.update(f_ringdown=m.freq, q_ringdown=m.q, r_peak_narrow=float(R[i]), f_peak_narrow=float(f[i]),
                       decay=bool(r.meta.get("decay_reached")), q_fwhm_narrow=float(Rz.fwhm_q(f, R, i)))
        res[f"{pitch:g}"] = rec
        print(pitch, {k: v for k, v in rec.items() if k != "harminv"}, flush=True)
        (REPO / "results" / "tier2" / "highq_summary.json").write_text(json.dumps(res, indent=1))
