"""Tier 2 primary study: how pitch and polarisation shape the R/T spectrum of a cylinder array.

Normal incidence, circular dielectric cylinders (R = 1, eps_r = 4: the scatterer verified in Tier 1)
on a one-dimensional lattice of pitch Lambda.  Results: ``results/tier2/pitch``.

    python experiments/tier2/pitch_study.py [--resolution 40]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load

REPO = Path(__file__).resolve().parents[2]
PITCHES = [2.2, 2.4, 2.6, 2.8, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]


def configs(resolution: float = 40.0):
    return [
        PeriodicConfig(
            name=f"pitch_{pol}_L{p:g}", polarization=pol, pitch=p, radius=1.0, eps_r=4.0, ky=0.0,
            resolution=resolution, f_min=0.05, f_max=0.60, n_freq=551, decay_tol=1e-7, max_run_time=1500.0,
        )
        for pol in ("TM", "TE")
        for p in PITCHES
    ]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolution", type=float, default=40.0)
    ap.add_argument("--out", type=Path, default=REPO / "results" / "tier2" / "pitch")
    args = ap.parse_args()
    for cfg in configs(args.resolution):
        r = run_periodic_or_load(cfg, str(args.out))
        print(f"{cfg.name:20s} wall {r.meta['wall_time_total_s']:6.1f}s decay={r.meta['decay_reached']}", flush=True)
