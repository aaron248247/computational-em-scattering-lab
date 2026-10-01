import numpy as np
import pytest
from em_rf_lab.config import ConfigError
from em_rf_lab.periodic import PeriodicConfig, diffraction_orders, incidence_angle, slab_rt, wood_anomaly_frequencies


@pytest.mark.parametrize("pol", ["TM", "TE"])
@pytest.mark.parametrize("ky", [0.0, 0.8, 1.7])
def test_slab_energy_conservation(pol, ky):
    f = np.linspace(0.3, 0.9, 50)
    R, T = slab_rt(f, 2.0, 1.3, ky, pol)
    np.testing.assert_allclose(R + T, 1.0, atol=1e-12)


def test_slab_normal_incidence_airy_closed_form():
    """T = 1/(1 + ((n-1/n)/2)^2 sin^2(2 pi f n d)) at normal incidence (independent formula)."""
    n, d = 2.0, 1.1
    f = np.linspace(0.1, 1.2, 80)
    _, T = slab_rt(f, n, d, 0.0, "TM")
    expect = 1.0 / (1.0 + ((n - 1 / n) / 2) ** 2 * np.sin(2 * np.pi * f * n * d) ** 2)
    np.testing.assert_allclose(T, expect, rtol=1e-12)
    _, Tp = slab_rt(f, n, d, 0.0, "TE")
    np.testing.assert_allclose(Tp, expect, rtol=1e-12)   # polarisation-degenerate at normal incidence


def test_slab_half_wave_and_quarter_wave():
    n, d = 2.0, 1.0
    f_half = 1.0 / (2 * n * d)            # delta = pi  -> transparent
    assert slab_rt(f_half, n, d)[1][0] == pytest.approx(1.0, abs=1e-12)
    f_q = 1.0 / (4 * n * d)               # delta = pi/2 -> reflectance maximum
    assert slab_rt(f_q, n, d)[0][0] == pytest.approx(((n**2 - 1) / (n**2 + 1)) ** 2, rel=1e-12)


def test_slab_brewster_p_pol():
    """p-polarised reflection from a single interface vanishes at Brewster; a slab is fully transmitting there."""
    n = 1.5
    theta_b = np.arctan(n)
    f = 0.37
    ky = 2 * np.pi * f * np.sin(theta_b)
    R, T = slab_rt(f, n, 0.9, ky, "TE")
    assert R[0] == pytest.approx(0.0, abs=1e-12) and T[0] == pytest.approx(1.0, abs=1e-12)
    Rs, _ = slab_rt(f, n, 0.9, ky, "TM")
    assert Rs[0] > 0.05


def test_slab_evanescent_is_nan():
    R, T = slab_rt(0.1, 2.0, 1.0, ky=2.0, polarization="TM")   # k0 = 0.63 < ky
    assert np.isnan(R[0]) and np.isnan(T[0])


def test_diffraction_orders_and_wood():
    assert diffraction_orders(0.3, 3.0) == [0]            # 1/pitch = 0.333 > 0.3
    assert diffraction_orders(0.4, 3.0) == [-1, 0, 1]
    w = wood_anomaly_frequencies(3.0, 0.0)
    assert w[0] == (-1, pytest.approx(1 / 3.0)) or w[0][1] == pytest.approx(1 / 3.0)
    # oblique: the -1 order emerges at f = (2pi/pitch - ky)/(2pi)
    w = dict(wood_anomaly_frequencies(3.0, 0.5))
    assert w[-1] == pytest.approx((2 * np.pi / 3.0 - 0.5) / (2 * np.pi))


def test_incidence_angle():
    th = incidence_angle(np.array([0.05, 0.2]), ky=0.5 * 2 * np.pi * 0.2)
    assert np.isnan(th[0]) and th[1] == pytest.approx(np.arcsin(0.5))


@pytest.mark.parametrize("bad", [dict(structure="sphere"), dict(polarization="X"), dict(pitch=1.5, radius=1.0),
                                 dict(f_min=0.5, f_max=0.2), dict(monitor_gap=6.0), dict(resolution=2), dict(n_freq=1)])
def test_periodic_config_validation(bad):
    with pytest.raises(ConfigError):
        PeriodicConfig(**bad)


def test_periodic_config_hash_and_roundtrip():
    c = PeriodicConfig()
    assert c.config_hash() == PeriodicConfig(name="other").config_hash()
    assert c.config_hash() != c.replace(pitch=3.1).config_hash()
    assert PeriodicConfig.from_dict(c.to_dict()) == c
    assert c.source_x < -c.monitor_x < -c.half_extent
