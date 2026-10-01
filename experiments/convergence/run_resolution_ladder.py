"""Resolution ladder: FDTD vs. exact Mie series as the grid is refined.

Usage
-----
    python experiments/convergence/run_resolution_ladder.py                 # full ladder
    python experiments/convergence/run_resolution_ladder.py --resolutions 10 15 20
    python experiments/convergence/run_resolution_ladder.py --quick         # small, for a smoke test

Results are cached by configuration hash under ``results/convergence/``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from em_rf_lab.config import SimulationConfig
from em_rf_lab.simulation import run_or_load

REPO = Path(__file__).resolve().parents[2]
DEFAULT_RESOLUTIONS = [10, 15, 20, 30, 40, 60, 80, 120]


def ladder_config(resolution: float, polarization: str = "TM", **overrides) -> SimulationConfig:
    return SimulationConfig(name=f"ladder_{polarization}_res{resolution:g}", resolution=resolution, polarization=polarization, **overrides)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolutions", type=float, nargs="+", default=DEFAULT_RESOLUTIONS)
    ap.add_argument("--polarization", choices=["TM", "TE"], default="TM")
    ap.add_argument("--out", type=Path, default=REPO / "results" / "convergence")
    ap.add_argument("--force", action="store_true", help="ignore cached results")
    ap.add_argument("--quick", action="store_true", help="coarse/cheap smoke-test settings")
    args = ap.parse_args()

    overrides = {}
    if args.quick:
        overrides = dict(n_freq=9, n_angles=31)
        args.resolutions = [r for r in args.resolutions if r <= 20]

    for res in args.resolutions:
        cfg = ladder_config(res, args.polarization, **overrides)
        r = run_or_load(cfg, str(args.out), force=args.force)
        print(
            f"res {res:>5g}  wall {r.meta['wall_time_total_s']:8.1f} s  "
            f"t_sim {r.meta['sim_time_scattering']:6.1f}  decay_ok={r.meta['decay_reached_scattering']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
