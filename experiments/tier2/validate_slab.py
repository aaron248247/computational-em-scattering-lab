"""Tier 2 verification: periodic-cell R/T workflow against the exact slab (Airy) solution.

A homogeneous slab has no pattern, so any deviation from the transfer-matrix result is a
property of the numerical workflow (Bloch boundaries, oblique source, flux normalisation,
PML), not of the physics.  Results: ``results/tier2/slab``.

    python experiments/tier2/validate_slab.py
"""

from __future__ import annotations

from pathlib import Path

from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load

REPO = Path(__file__).resolve().parents[2]
OUT = str(REPO / "results" / "tier2" / "slab")
BASE = dict(structure="slab", thickness=1.0, eps_r=4.0, pitch=1.0, f_min=0.3, f_max=0.9, n_freq=61, decay_tol=1e-9)


def configs():
    cfgs = []
    for pol in ("TM", "TE"):
        for res in (20, 30, 40, 60, 80):                      # normal-incidence resolution ladder
            cfgs.append(PeriodicConfig(name=f"slab_{pol}_normal_res{res}", polarization=pol, ky=0.0, resolution=res, **BASE))
        for ky in (0.5, 1.0, 1.5):                            # oblique incidence at fixed ky
            cfgs.append(PeriodicConfig(name=f"slab_{pol}_ky{ky:g}_res40", polarization=pol, ky=ky, resolution=40,
                                       max_run_time=300.0, **BASE))
    return cfgs


if __name__ == "__main__":
    for cfg in configs():
        r = run_periodic_or_load(cfg, OUT)
        print(f"{cfg.name:30s} wall {r.meta['wall_time_total_s']:6.1f}s", flush=True)
