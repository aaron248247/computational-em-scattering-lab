"""Resonance characterisation from reflectance spectra (pure NumPy/SciPy).

A single isolated resonance coupled to one radiation channel gives the
coupled-mode (Fano) form for the complex reflection amplitude::

    r(f) = r_b + A / (f - f0 + i*gamma)

where ``r_b`` is a slowly-varying background, ``A`` a complex coupling and
``Q = f0 / (2*gamma)`` the radiative quality factor (``gamma`` = half-width at
half maximum of the Lorentzian pole).  ``R(f) = |r(f)|^2``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

__all__ = ["ResonanceFit", "cmt_reflectance", "fit_resonance", "find_reflection_peaks", "fwhm_q"]


def cmt_reflectance(f, f0, gamma, rb, a):
    """``|rb + a/(f - f0 + i gamma)|^2`` for complex ``rb``, ``a``."""
    f = np.asarray(f, dtype=float)
    return np.abs(rb + a / (f - f0 + 1j * gamma)) ** 2


@dataclass(frozen=True)
class ResonanceFit:
    f0: float
    gamma: float
    rb: complex
    a: complex
    rms_residual: float
    n_points: int

    @property
    def q(self) -> float:
        return self.f0 / (2.0 * self.gamma)

    @property
    def fwhm(self) -> float:
        return 2.0 * self.gamma


def find_reflection_peaks(f, R, min_height=0.3, min_prominence=0.15):
    """Indices of reflectance maxima."""
    idx, _ = find_peaks(np.asarray(R), height=min_height, prominence=min_prominence)
    return idx


def fwhm_q(f, R, i_peak: int, baseline: float = 0.0) -> float:
    """Model-free ``Q = f_peak / FWHM`` from linear interpolation of half-maximum crossings.

    Half maximum is taken between ``baseline`` (default 0, a dark background) and the peak.
    """
    f = np.asarray(f, dtype=float)
    R = np.asarray(R, dtype=float)
    half = 0.5 * (R[i_peak] + baseline)
    lo = i_peak
    while lo > 0 and R[lo] > half:
        lo -= 1
    hi = i_peak
    while hi < len(R) - 1 and R[hi] > half:
        hi += 1
    if lo == 0 or hi == len(R) - 1:
        return float("nan")
    f_lo = f[lo] + (half - R[lo]) * (f[lo + 1] - f[lo]) / (R[lo + 1] - R[lo])
    f_hi = f[hi - 1] + (half - R[hi - 1]) * (f[hi] - f[hi - 1]) / (R[hi] - R[hi - 1])
    return float(f[i_peak] / (f_hi - f_lo))


def fit_resonance(f, R, f_guess: float, half_window: float, gamma_guess: float | None = None) -> ResonanceFit:
    """Least-squares CMT fit of ``R(f)`` in ``[f_guess - w, f_guess + w]``."""
    f = np.asarray(f, dtype=float)
    R = np.asarray(R, dtype=float)
    sel = (f >= f_guess - half_window) & (f <= f_guess + half_window)
    if sel.sum() < 8:
        raise ValueError("fit window contains fewer than 8 samples")
    fs, Rs = f[sel], R[sel]
    g0 = gamma_guess or max(half_window / 6.0, 1e-4)

    def model(p):
        f0, g, rbr, rbi, ar, ai = p
        return cmt_reflectance(fs, f0, g, rbr + 1j * rbi, ar + 1j * ai)

    best = None
    for scale in (0.5, 1.0, 2.0):
        for ph in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2):
            amp = g0 * scale
            p0 = [f_guess, g0 * scale, 0.0, 0.0, amp * np.cos(ph), amp * np.sin(ph)]
            try:
                sol = least_squares(lambda p: model(p) - Rs, p0, bounds=([f_guess - half_window, 1e-6, -2, -2, -1, -1], [f_guess + half_window, half_window, 2, 2, 1, 1]))
            except Exception:
                continue
            if best is None or sol.cost < best.cost:
                best = sol
    if best is None:
        raise RuntimeError("resonance fit failed")
    f0, g, rbr, rbi, ar, ai = best.x
    rms = float(np.sqrt(np.mean((model(best.x) - Rs) ** 2)))
    return ResonanceFit(float(f0), float(g), complex(rbr, rbi), complex(ar, ai), rms, int(sel.sum()))
