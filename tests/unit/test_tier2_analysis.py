import numpy as np
import pytest
from em_rf_lab import resonance as Rz
from em_rf_lab.tier2_analysis import band_edges, energy_residual, summarise_peaks


def _lorentz(f, f0, q, peak=1.0):
    g = f0 / (2 * q)
    return peak * g**2 / ((f - f0) ** 2 + g**2)


def test_band_edges_of_lorentzian():
    f0, q = 0.3, 40.0
    f = np.linspace(0.2, 0.4, 4001)
    R = _lorentz(f, f0, q)
    e = band_edges(f, R, int(np.argmax(R)), threshold=0.5)
    assert e.width == pytest.approx(f0 / q, rel=2e-3)             # FWHM
    e9 = band_edges(f, R, int(np.argmax(R)), threshold=0.9)
    g = f0 / (2 * q)
    assert e9.width == pytest.approx(2 * g * np.sqrt(1 / 0.9 - 1), rel=5e-3)
    assert band_edges(f, 0.5 * R, int(np.argmax(R)), 0.9) is None


def test_summarise_recovers_q_and_flags():
    f = np.linspace(0.1, 0.5, 2001)
    R = _lorentz(f, 0.25, 30.0) + _lorentz(f, 0.4, 80.0, 0.8)
    T = 1.0 - R
    peaks = summarise_peaks(f, R, T)
    assert [round(p.f_peak, 3) for p in peaks] == [0.25, 0.4]
    assert peaks[0].q_fwhm == pytest.approx(30.0, rel=0.03) and peaks[1].q_fit == pytest.approx(80.0, rel=0.05)
    assert all("energy_not_conserved" not in p.flags for p in peaks)
    bad = summarise_peaks(f, R, T - 0.2)
    assert all("energy_not_conserved" in p.flags for p in bad)


def test_undersampled_flag():
    f = np.linspace(0.1, 0.5, 61)            # 0.0067 spacing vs FWHM 0.0025
    R = _lorentz(f, 0.3, 120.0)
    peaks = summarise_peaks(f, R, 1 - R, min_prominence=0.05)
    assert peaks and "undersampled_line" in peaks[0].flags


def test_energy_residual():
    np.testing.assert_allclose(energy_residual([0.3], [0.69]), [0.01])
