"""Post-processing for periodic-layer (FSS / metasurface) spectra."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import resonance as Rz

__all__ = ["BandEdges", "PeakSummary", "band_edges", "summarise_peaks", "energy_residual"]


@dataclass(frozen=True)
class BandEdges:
    f_lo: float
    f_hi: float

    @property
    def width(self) -> float:
        return self.f_hi - self.f_lo

    @property
    def fractional(self) -> float:
        return self.width / (0.5 * (self.f_hi + self.f_lo))


def band_edges(f, R, i_peak: int, threshold: float = 0.9):
    """Contiguous band around ``i_peak`` where ``R >= threshold`` (linear interpolation at the edges).

    Returns ``None`` if the peak itself is below ``threshold``.
    """
    f = np.asarray(f, dtype=float)
    R = np.asarray(R, dtype=float)
    if R[i_peak] < threshold:
        return None
    lo = i_peak
    while lo > 0 and R[lo] >= threshold:
        lo -= 1
    hi = i_peak
    while hi < len(R) - 1 and R[hi] >= threshold:
        hi += 1
    f_lo = f[lo] if R[lo] >= threshold else f[lo] + (threshold - R[lo]) * (f[lo + 1] - f[lo]) / (R[lo + 1] - R[lo])
    f_hi = f[hi] if R[hi] >= threshold else f[hi - 1] + (threshold - R[hi - 1]) * (f[hi] - f[hi - 1]) / (R[hi] - R[hi - 1])
    return BandEdges(float(f_lo), float(f_hi))


@dataclass
class PeakSummary:
    f_peak: float
    r_peak: float
    q_fwhm: float                  # model-free
    q_fit: float                   # coupled-mode fit (nan if the fit failed)
    f_fit: float
    fit_rms: float
    band90: BandEdges | None
    max_abs_A_near_peak: float     # energy-conservation residual within the peak window
    flags: list[str] = field(default_factory=list)


def energy_residual(R, T) -> np.ndarray:
    """``|1 - R - T|`` (zero for a lossless, converged simulation)."""
    return np.abs(1.0 - np.asarray(R) - np.asarray(T))


def summarise_peaks(f, R, T, min_height=0.3, min_prominence=0.1, a_tol=0.01) -> list[PeakSummary]:
    """Characterise every reflection peak: frequency, height, Q (two ways), >=90% band, conservation."""
    f = np.asarray(f, dtype=float)
    R = np.asarray(R, dtype=float)
    T = np.asarray(T, dtype=float)
    out: list[PeakSummary] = []
    df = float(f[1] - f[0])
    for i in Rz.find_reflection_peaks(f, R, min_height, min_prominence):
        q_fwhm = Rz.fwhm_q(f, R, int(i))
        half = max(6 * df, (f[i] / q_fwhm) * 1.5 if np.isfinite(q_fwhm) and q_fwhm > 0 else 12 * df)
        flags: list[str] = []
        try:
            fit = Rz.fit_resonance(f, R, f_guess=float(f[i]), half_window=half, gamma_guess=half / 6)
            qf, ff, rms = fit.q, fit.f0, fit.rms_residual
            if rms > 0.02:
                flags.append("poor_fit")
        except Exception:
            qf, ff, rms = float("nan"), float("nan"), float("nan")
            flags.append("fit_failed")
        window = np.abs(f - f[i]) <= half
        a_res = float(np.max(energy_residual(R[window], T[window])))
        if a_res > a_tol:
            flags.append("energy_not_conserved")
        if np.isfinite(q_fwhm) and (f[1] - f[0]) > 0.4 * f[i] / q_fwhm:
            flags.append("undersampled_line")
        out.append(PeakSummary(float(f[i]), float(R[i]), float(q_fwhm), float(qf), float(ff), float(rms), band_edges(f, R, int(i)), a_res, flags))
    return out
