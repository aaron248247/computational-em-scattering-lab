"""Unit tests for the exact cylinder-scattering reference (em_rf_lab.analytical)."""

import numpy as np
import pytest
from scipy import special

from em_rf_lab import analytical as A

POLS = ["TM", "TE"]
MATS = [("dielectric", 4.0), ("dielectric", 2.25), ("pec", 1.0)]


@pytest.mark.parametrize("pol", POLS)
@pytest.mark.parametrize("material,eps", MATS)
@pytest.mark.parametrize("x", [0.3, 1.0, 2.7, 6.1, 11.0])
def test_optical_theorem_lossless(pol, material, eps, x):
    """Lossless scatterer: extinction width equals scattering width (energy conservation)."""
    s = A.total_scattering_width(x, 1.0, eps, pol, material)
    e = A.total_extinction_width(x, 1.0, eps, pol, material)
    assert s > 0
    assert e == pytest.approx(s, rel=1e-12)


@pytest.mark.parametrize("pol", POLS)
def test_lossy_cylinder_absorbs(pol):
    """Complex permittivity: extinction exceeds scattering (positive absorption)."""
    x = 2.0
    s = A.total_scattering_width(x, 1.0, 4.0 + 0.5j, pol)
    e = A.total_extinction_width(x, 1.0, 4.0 + 0.5j, pol)
    assert e > s * (1 + 1e-3)


@pytest.mark.parametrize("pol", POLS)
@pytest.mark.parametrize("material,eps", MATS)
def test_echo_width_integrates_to_total_width(pol, material, eps):
    """(1/2pi) * integral of W over the circle equals sigma_s (independent quadrature)."""
    x = 3.3
    phi = np.linspace(0, 2 * np.pi, 4001)[:-1]
    w = A.echo_width(phi, x, 1.0, eps, pol, material)
    assert w.mean() == pytest.approx(A.total_scattering_width(x, 1.0, eps, pol, material), rel=1e-10)


@pytest.mark.parametrize("pol", POLS)
def test_echo_width_symmetry_about_axis(pol):
    phi = np.linspace(0.1, 3.0, 17)
    np.testing.assert_allclose(A.echo_width(phi, 2.2, 1.0, 4.0, pol), A.echo_width(-phi, 2.2, 1.0, 4.0, pol), rtol=1e-12)


@pytest.mark.parametrize("pol", POLS)
@pytest.mark.parametrize("x", [0.8, 2.5, 5.0])
def test_coefficients_match_direct_boundary_solve(pol, x):
    """Closed-form a_n agrees with a direct numerical solve of the 2x2 matching system."""
    eps = 4.0
    m = np.sqrt(eps)
    n_max = 25
    a = A.mie_coefficients(x, eps, pol, n_max=n_max)
    for n in (0, 1, 2, 5, 9):
        # unknowns (a_n, c_n):  J + a H = c J(mx);  J' + a H' = g c J'(mx)   [g=m (TM), 1/m (TE)]
        g = m if pol == "TM" else 1.0 / m
        M = np.array(
            [
                [special.hankel1(n, x), -special.jv(n, m * x)],
                [special.h1vp(n, x), -g * special.jvp(n, m * x)],
            ]
        )
        rhs = -np.array([special.jv(n, x), special.jvp(n, x)])
        a_n, _ = np.linalg.solve(M, rhs)
        assert a[n] == pytest.approx(a_n, rel=1e-10, abs=1e-14)


@pytest.mark.parametrize("pol", POLS)
def test_near_field_satisfies_interface_conditions(pol):
    """F continuous across rho=R; (1/eps)^(TE) dF/drho continuous. Checked by finite differences."""
    x, eps, h = 3.0, 4.0, 1e-5
    th = np.linspace(0.0, np.pi, 9)

    def field(r):
        return A.near_field(r * np.cos(th), r * np.sin(th), x, 1.0, eps, pol)

    f_in, f_out = field(1 - h), field(1 + h)
    np.testing.assert_allclose(f_in, f_out, atol=5e-4)
    d_in = (field(1 - h) - field(1 - 3 * h)) / (2 * h)
    d_out = (field(1 + 3 * h) - field(1 + h)) / (2 * h)
    weight_in = 1.0 / eps if pol == "TE" else 1.0
    np.testing.assert_allclose(weight_in * d_in, d_out, atol=2e-3)


def test_pec_field_vanishes_on_surface_tm():
    th = np.linspace(0, np.pi, 7)
    f = A.near_field(np.cos(th), np.sin(th), 2.0, 1.0, polarization="TM", material="pec")
    np.testing.assert_allclose(np.abs(f), 0.0, atol=1e-10)


def test_rayleigh_limit_tm_dielectric():
    """Small-cylinder TM limit: sigma_s -> (pi^2/4)(eps-1)^2 k^3 R^4."""
    x, eps = 0.01, 4.0
    expected = (np.pi**2 / 4) * (eps - 1) ** 2 * x**3
    assert A.total_scattering_width(x, 1.0, eps, "TM") == pytest.approx(expected, rel=2e-3)


def test_rayleigh_scaling_exponent():
    """log-log slope of sigma_s vs ka is 3 for a small dielectric cylinder (TM)."""
    s1 = A.total_scattering_width(0.01, 1.0, 4.0, "TM")
    s2 = A.total_scattering_width(0.02, 1.0, 4.0, "TM")
    assert np.log2(s2 / s1) == pytest.approx(3.0, abs=0.02)


def test_pec_geometric_optics_limit():
    """Large PEC cylinder: Q_sca = sigma_s/(2R) -> 2 (extinction paradox)."""
    for pol in POLS:
        q = A.scattering_efficiency(80.0, polarization=pol, material="pec")
        assert q == pytest.approx(2.0, abs=0.12)


def test_dielectric_reduces_to_vacuum_for_eps_one():
    assert A.total_scattering_width(3.0, 1.0, 1.0, "TM") == pytest.approx(0.0, abs=1e-20)


def test_series_truncation_converged():
    """Adding orders beyond the default changes nothing measurable."""
    x = 7.5
    n_default = A.default_n_max(x, 2.0)
    w0 = A.echo_width(1.0, x, 1.0, 4.0, "TM", n_max=n_default)
    w1 = A.echo_width(1.0, x, 1.0, 4.0, "TM", n_max=n_default + 30)
    assert w0 == pytest.approx(w1, rel=1e-13)


def test_radius_scaling():
    """W scales linearly with length when x = kR is held fixed."""
    w1 = A.echo_width(1.0, 3.0, radius=1.0)
    w2 = A.echo_width(1.0, 3.0, radius=2.5)
    assert w2 == pytest.approx(2.5 * w1, rel=1e-12)


def test_backscatter_matches_echo_width_at_pi():
    assert A.backscatter_width(2.0) == pytest.approx(A.echo_width(np.pi, 2.0))


def test_regression_reference_values():
    """Pinned values (computed and cross-checked via the optical theorem and direct 2x2 solve)."""
    assert A.total_scattering_width(3.0, 1.0, 4.0, "TM") == pytest.approx(5.165944073648, rel=1e-10)
    assert A.total_scattering_width(3.0, 1.0, 4.0, "TE") == pytest.approx(3.391626, rel=1e-6)
    assert A.total_scattering_width(1.0, 1.0, material="pec", polarization="TM") == pytest.approx(5.913114, rel=1e-6)


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_invalid_size_parameter(bad):
    with pytest.raises(ValueError):
        A.mie_coefficients(bad)


def test_invalid_polarization_and_material():
    with pytest.raises(ValueError):
        A.mie_coefficients(1.0, polarization="TEM")
    with pytest.raises(ValueError):
        A.mie_coefficients(1.0, material="plasma")
