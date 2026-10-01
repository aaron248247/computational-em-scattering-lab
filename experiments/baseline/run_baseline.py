"""Baseline benchmark run (TM, eps_r = 4, res 40) with a stored near-field map at ka = 3.

    python experiments/baseline/run_baseline.py
"""

from __future__ import annotations

from pathlib import Path

from em_rf_lab.config import SimulationConfig
from em_rf_lab.simulation import run_or_load

REPO = Path(__file__).resolve().parents[2]

if __name__ == "__main__":
    cfg = SimulationConfig(name="baseline_TM_res40", resolution=40.0)
    r = run_or_load(cfg, str(REPO / "results" / "baseline"), map_ka=3.0)
    print(f"baseline done: wall {r.meta['wall_time_total_s']:.1f} s")
