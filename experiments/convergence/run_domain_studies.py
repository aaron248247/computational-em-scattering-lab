"""One-factor-at-a-time numerical-setup studies at fixed resolution.

Varies, one at a time around the production configuration: PML thickness,
Gaussian beam waist (finite-beam effect), near-to-far box size, run duration
(time after the source) and source bandwidth.  Results: ``results/domain/``.

    python experiments/convergence/run_domain_studies.py [--resolution 30]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from em_rf_lab.config import SimulationConfig
from em_rf_lab.simulation import run_or_load

REPO = Path(__file__).resolve().parents[2]


def study_configs(resolution: float) -> dict[str, list[SimulationConfig]]:
    def base(**kw):
        name = kw.pop("name")
        return SimulationConfig(name=name, resolution=resolution, **kw)

    studies: dict[str, list[SimulationConfig]] = {"pml": [], "beam_waist": [], "n2f_half_width": [], "run_time": [], "source_bw": []}
    for pml in (2.0, 3.0, 5.0, 8.0):
        studies["pml"].append(base(name=f"dom_pml{pml:g}", pml=pml, cell_x=2 * (pml + 5.0)))
    for w in (5.0, 10.0, 20.0, 30.0):
        studies["beam_waist"].append(base(name=f"dom_waist{w:g}", beam_waist=w, cell_y=max(56.0, 2 * (5.0 + 1.25 * w))))
    for h in (2.0, 2.5, 3.5, 5.0):
        studies["n2f_half_width"].append(base(name=f"dom_n2f{h:g}", n2f_half_width=h, cell_x=max(20.0, 2 * (5.0 + h + 2.0))))
    for T in (50.0, 110.0, 250.0, 500.0, 1000.0):
        # decay_tol = 1 is satisfied immediately, so the run lasts exactly min_run_time after the source
        studies["run_time"].append(base(name=f"dom_T{T:g}", decay_tol=1.0, min_run_time=T, max_run_time=T + 1.0))
    for f in (0.75, 1.5, 3.0):
        studies["source_bw"].append(base(name=f"dom_bw{f:g}", source_fwidth_factor=f))
    return studies


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolution", type=float, default=30.0)
    ap.add_argument("--out", type=Path, default=REPO / "results" / "domain")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    for study, cfgs in study_configs(args.resolution).items():
        for cfg in cfgs:
            r = run_or_load(cfg, str(args.out), force=args.force)
            print(f"{study:15s} {cfg.name:18s} wall {r.meta['wall_time_total_s']:6.1f}s", flush=True)


if __name__ == "__main__":
    main()
