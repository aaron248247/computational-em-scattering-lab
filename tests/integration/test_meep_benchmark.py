"""FDTD-vs-Mie integration tests (Tier 1).  Cheap ones run by default; `-m slow` adds heavier ones."""

import glob
from pathlib import Path

import numpy as np
import pytest

from em_rf_lab import convergence as C
from em_rf_lab.config import SimulationConfig
from em_rf_lab.io import load_run
from em_rf_lab.simulation import run_scattering
from em_rf_lab.validation import compare_to_mie

REPO = Path(__file__).resolve().parents[2]
SMALL = dict(n_freq=5, n_angles=31, ka_min=2.0, ka_max=6.0)


@pytest.mark.meep
def test_small_run_is_consistent_and_close_to_mie():
    r = run_scattering(SimulationConfig(name="t", resolution=15, **SMALL))
    assert r.meta["decay_reached_scattering"] and r.meta["decay_reached_reference"]
    assert np.all(np.abs(r["flux_norm_check"] - 1.0) < 0.15)            # flux calibration (factor 2 convention)
    # two independent estimators of sigma_s agree (near-to-far integral vs scattered flux)
    np.testing.assert_allclose(r["sigma_s_flux"], r["sigma_s_n2f"], rtol=0.06)
    c = compare_to_mie(r)
    assert c.rms_rel_l2_angle < 0.30                                    # coarse grid: loose bound
    assert np.all(c.sigma_ref > 0)


@pytest.mark.meep
def test_error_decreases_with_resolution():
    e = [compare_to_mie(run_scattering(SimulationConfig(name="t", resolution=res, **SMALL))).rms_rel_l2_angle for res in (12, 24)]
    assert e[1] < 0.6 * e[0]


@pytest.mark.meep
@pytest.mark.parametrize("pol", ["TM", "TE"])
def test_mirror_symmetry_reproduces_full_domain(pol):
    kw = dict(name="t", resolution=10, polarization=pol, n_freq=3, n_angles=11, cell_y=30.0, beam_waist=10.0)
    a = run_scattering(SimulationConfig(use_symmetry=False, **kw))
    b = run_scattering(SimulationConfig(use_symmetry=True, **kw))
    np.testing.assert_allclose(a["echo_width"], b["echo_width"], rtol=1e-8)
    np.testing.assert_allclose(a["sigma_s_flux"], b["sigma_s_flux"], rtol=1e-8)


@pytest.mark.meep
def test_te_polarisation_runs_and_matches_te_mie():
    r = run_scattering(SimulationConfig(name="t", resolution=15, polarization="TE", **SMALL))
    c = compare_to_mie(r)
    assert c.rms_rel_l2_angle < 0.35


@pytest.mark.results
def test_committed_ladder_results_meet_validation_criteria():
    paths = glob.glob(str(REPO / "results" / "convergence" / "ladder_TM_*"))
    if len(paths) < 4:
        pytest.skip("ladder results not present")
    runs = sorted((load_run(p) for p in paths), key=lambda r: r.config.resolution)
    res = np.array([r.config.resolution for r in runs])
    rms = np.array([compare_to_mie(r).rms_rel_l2_angle for r in runs])
    assert np.all(np.diff(rms) < 0), "error must decrease monotonically with resolution"
    fit = C.fit_order(res[res >= 20], rms[res >= 20])
    assert 1.5 < fit.p < 2.5, f"observed order {fit.p:.2f} outside the second-order band"
    assert rms[res == 40][0] < 0.05
    assert all(r.meta["decay_reached_scattering"] for r in runs)
