"""Error metrics for comparing numerical and reference curves.

Conventions
-----------
For a numerical array ``q`` and a reference ``ref`` (same shape):

* absolute error           ``|q - ref|``
* relative L2 error        ``||q - ref||_2 / ||ref||_2`` (optionally with
  trapezoidal weights ``x`` so the norm approximates the continuous L2 norm)
* maximum relative error   ``max |q - ref| / max |ref|``  (normalised by the peak
  of the reference so near-zeros of the reference do not dominate)
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "abs_error",
    "rel_error",
    "rel_l2_error",
    "max_abs_deviation",
    "max_rel_deviation",
    "successive_difference",
    "log10_ratio_db",
]


def _as_pair(q, ref):
    q = np.asarray(q)
    ref = np.asarray(ref)
    if q.shape != ref.shape:
        raise ValueError(f"shape mismatch: {q.shape} vs {ref.shape}")
    return q, ref


def abs_error(q, ref) -> np.ndarray:
    q, ref = _as_pair(q, ref)
    return np.abs(q - ref)


def rel_error(q, ref) -> np.ndarray:
    """Pointwise ``|q - ref| / |ref|``; raises if any reference value is zero."""
    q, ref = _as_pair(q, ref)
    if np.any(ref == 0):
        raise ValueError("reference contains zeros; pointwise relative error undefined")
    return np.abs(q - ref) / np.abs(ref)


def rel_l2_error(q, ref, x=None) -> float:
    """Relative L2 error, optionally integrating over abscissa ``x`` (trapezoid rule)."""
    q, ref = _as_pair(q, ref)
    d2 = np.abs(q - ref) ** 2
    r2 = np.abs(ref) ** 2
    if x is None:
        num, den = d2.sum(), r2.sum()
    else:
        x = np.asarray(x, dtype=float)
        if x.shape != q.shape:
            raise ValueError("x must have the same shape as the data")
        num, den = np.trapezoid(d2, x), np.trapezoid(r2, x)
    if den == 0:
        raise ValueError("reference has zero norm")
    return float(np.sqrt(num / den))


def max_abs_deviation(q, ref) -> float:
    return float(np.max(abs_error(q, ref)))


def max_rel_deviation(q, ref) -> float:
    """``max|q - ref| / max|ref|`` (peak-normalised)."""
    q, ref = _as_pair(q, ref)
    peak = np.max(np.abs(ref))
    if peak == 0:
        raise ValueError("reference is identically zero")
    return float(np.max(np.abs(q - ref)) / peak)


def successive_difference(q_h, q_h2) -> float:
    """Grid-convergence diagnostic ``eps_h = ||Q_h - Q_{h/2}|| / ||Q_{h/2}||`` (L2 over samples)."""
    return rel_l2_error(q_h, q_h2)


def log10_ratio_db(q, ref) -> np.ndarray:
    """``10 log10(q / ref)`` in dB for positive arrays (for intensity-ratio plots)."""
    q, ref = _as_pair(q, ref)
    if np.any(q <= 0) or np.any(ref <= 0):
        raise ValueError("arguments must be strictly positive")
    return 10.0 * np.log10(q / ref)
