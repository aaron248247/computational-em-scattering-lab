import numpy as np
import pytest
from em_rf_lab import metrics as M, convergence as C
from em_rf_lab.config import SimulationConfig
from em_rf_lab.farfield import echo_width_from_farfield, total_width_from_echo
from em_rf_lab.geometry import angle_grid
from em_rf_lab import analytical as A
from em_rf_lab.io import load_run, run_dir_name, save_run


def test_metrics_basic():
    ref = np.array([1.0, 2.0, 4.0])
    q = ref * 1.1
    assert M.rel_l2_error(q, ref) == pytest.approx(0.1)
    assert M.max_abs_deviation(q, ref) == pytest.approx(0.4)
    assert M.max_rel_deviation(q, ref) == pytest.approx(0.1)
    np.testing.assert_allclose(M.rel_error(q, ref), 0.1)
    assert M.rel_l2_error(ref, ref) == 0.0


def test_metrics_weighted_l2_and_errors():
    x = np.linspace(0, 1, 101)
    assert M.rel_l2_error(1.05 * np.ones_like(x), np.ones_like(x), x) == pytest.approx(0.05)
    with pytest.raises(ValueError):
        M.rel_l2_error([1, 2], [1, 2, 3])
    with pytest.raises(ValueError):
        M.rel_error([1.0], [0.0])
    with pytest.raises(ValueError):
        M.rel_l2_error([0.0], [0.0])
    assert M.log10_ratio_db([10.0], [1.0])[0] == pytest.approx(10.0)


def test_fit_order_recovers_exponent():
    res = np.array([10, 20, 40, 80.0])
    err = 3.0 * (1 / res) ** 2
    fit = C.fit_order(res, err)
    assert fit.p == pytest.approx(2.0) and fit.C == pytest.approx(3.0) and fit.r2 == pytest.approx(1.0)
    np.testing.assert_allclose(C.local_orders(res, err), 2.0)
    with pytest.raises(ValueError):
        C.fit_order(res[:1], err[:1])


def test_richardson_order():
    assert C.richardson_order(4.0, 1.0, 2.0) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        C.richardson_order(1.0, 0.0)


def test_echo_width_normalisation_roundtrip():
    rho, w = 1e5, 3.7
    f_inc = np.array([0.5 + 0.1j])
    f_s = np.sqrt(w / (2 * np.pi * rho)) * np.abs(f_inc)[0] * np.ones((1, 5)) * np.exp(1j * 0.3)
    np.testing.assert_allclose(echo_width_from_farfield(f_s, f_inc, rho), w)
    with pytest.raises(ValueError):
        echo_width_from_farfield(f_s, f_inc, -1.0)


def test_total_width_from_echo_matches_mie_sigma():
    ang = angle_grid(91)
    w = A.echo_width(ang, 3.0, 1.0, 4.0, "TM")[None, :]
    s = total_width_from_echo(w, ang)[0]
    assert s == pytest.approx(A.total_scattering_width(3.0, 1.0, 4.0, "TM"), rel=1e-9)


def test_save_load_roundtrip_and_missing(tmp_path):
    cfg = SimulationConfig(name="io test")
    arrays = {"ka": np.linspace(2, 6, 5), "z": np.array([1 + 2j])}
    p = save_run(tmp_path, cfg, arrays, {"note": "x"})
    assert p.name == run_dir_name(cfg) and " " not in p.name
    r = load_run(p)
    assert r.config == cfg and r.meta["note"] == "x" and "numpy" in r.meta
    np.testing.assert_array_equal(r["z"], arrays["z"])
    (p / "arrays.npz").unlink()
    with pytest.raises(FileNotFoundError):
        load_run(p)
