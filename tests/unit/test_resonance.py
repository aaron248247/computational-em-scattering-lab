import numpy as np
import pytest
from em_rf_lab import resonance as Rz


def _synthetic(f0, q, rb, a_phase, n=801, width=0.05):
    g = f0 / (2 * q)
    f = np.linspace(f0 - width, f0 + width, n)
    # choose |A| so that the lossless-like peak R reaches ~1 for rb = 0; generic otherwise
    a = g * np.exp(1j * a_phase)
    return f, Rz.cmt_reflectance(f, f0, g, rb, a), g


@pytest.mark.parametrize("q", [15.0, 40.0, 120.0])
@pytest.mark.parametrize("rb,ph", [(0.0, 0.0), (0.3 + 0.1j, 1.0), (-0.4j, 2.5)])
def test_fit_recovers_known_q(q, rb, ph):
    f0 = 0.3
    f, R, g = _synthetic(f0, q, rb, ph, width=max(0.05, 6 * f0 / (2 * q)))
    fit = Rz.fit_resonance(f, R, f_guess=f0 + 0.2 * g, half_window=max(0.03, 5 * g), gamma_guess=1.5 * g)
    assert fit.f0 == pytest.approx(f0, abs=2e-3 * f0)
    assert fit.q == pytest.approx(q, rel=0.03)
    assert fit.rms_residual < 1e-6


def test_fwhm_q_of_lorentzian():
    f0, q = 0.25, 30.0
    g = f0 / (2 * q)
    f = np.linspace(0.2, 0.3, 4001)
    R = g**2 / ((f - f0) ** 2 + g**2)
    i = int(np.argmax(R))
    assert Rz.fwhm_q(f, R, i) == pytest.approx(q, rel=2e-3)


def test_peak_finder_and_window_error():
    f = np.linspace(0, 1, 500)
    R = np.exp(-((f - 0.4) / 0.01) ** 2) + 0.5 * np.exp(-((f - 0.7) / 0.02) ** 2)
    assert list(f[Rz.find_reflection_peaks(f, R)].round(2)) == [0.4, 0.7]
    with pytest.raises(ValueError):
        Rz.fit_resonance(f, R, 0.4, 1e-4)
