"""Checks on the committed Tier 2 result files (no MEEP needed; run with ``pytest -m results``)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.results


def _load(name):
    p = REPO / "results" / "tier2" / name
    if not p.exists():
        pytest.skip(f"{name} not present")
    return json.loads(p.read_text())


def test_resonance_metrics_converged_across_resolution_and_domain():
    s = _load("convergence_summary.json")
    for case in s.values():
        for study in case.values():
            f = [v["f_fit"] for v in study.values()]
            q = [v["q_fit"] for v in study.values()]
            assert max(f) / min(f) - 1 < 1e-3
            assert max(q) / min(q) - 1 < 5e-3


def test_lossless_te_band_fully_reflecting():
    d = _load("lowest_resonance_summary.json")
    for key, v in d.items():
        if key.startswith("TE"):
            assert v["r_peak"] > 0.99 and v["max_abs_A"] < 0.01


def test_harminv_matches_line_shape_q_for_tm_pitch3():
    d = _load("lowest_resonance_summary.json")["TM_3"]
    h = [m for m in d["harminv"] if abs(m["freq"] - d["f_fit"]) < 0.005]
    assert h and abs(h[0]["q"] / d["q_fit"] - 1) < 0.15
