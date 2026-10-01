"""Convergence-rate analysis for grid-refinement studies."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["OrderFit", "fit_order", "local_orders", "richardson_order"]


@dataclass(frozen=True)
class OrderFit:
    """Least-squares fit ``err ~ C * h^p`` in log-log space."""

    p: float
    C: float
    r2: float
    n_points: int

    def predict(self, h):
        return self.C * np.asarray(h, dtype=float) ** self.p


def fit_order(resolution, err) -> OrderFit:
    """Fit ``err = C h^p`` with ``h = 1/resolution``; requires >= 2 positive points."""
    res = np.asarray(resolution, dtype=float)
    err = np.asarray(err, dtype=float)
    if res.shape != err.shape or res.ndim != 1:
        raise ValueError("resolution and err must be 1-D arrays of equal length")
    mask = (err > 0) & np.isfinite(err) & (res > 0)
    if mask.sum() < 2:
        raise ValueError("need at least two positive finite error values")
    lh = np.log(1.0 / res[mask])
    le = np.log(err[mask])
    p, logC = np.polyfit(lh, le, 1)
    pred = p * lh + logC
    ss_res = float(np.sum((le - pred) ** 2))
    ss_tot = float(np.sum((le - le.mean()) ** 2))
    r2 = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return OrderFit(p=float(p), C=float(np.exp(logC)), r2=r2, n_points=int(mask.sum()))


def local_orders(resolution, err) -> np.ndarray:
    """Pairwise observed orders ``log(e_i/e_{i+1}) / log(res_{i+1}/res_i)``."""
    res = np.asarray(resolution, dtype=float)
    err = np.asarray(err, dtype=float)
    if res.shape != err.shape or res.size < 2:
        raise ValueError("need >= 2 matching points")
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(err[:-1] / err[1:]) / np.log(res[1:] / res[:-1])


def richardson_order(d_coarse: float, d_fine: float, ratio: float = 2.0) -> float:
    """Observed order from successive differences on three grids.

    With ``d_coarse = ||Q_h - Q_{h/r}||`` and ``d_fine = ||Q_{h/r} - Q_{h/r^2}||``
    the asymptotic order is ``p = log(d_coarse / d_fine) / log(r)``.
    """
    if d_coarse <= 0 or d_fine <= 0:
        raise ValueError("differences must be positive")
    if ratio <= 1:
        raise ValueError("refinement ratio must exceed 1")
    return float(np.log(d_coarse / d_fine) / np.log(ratio))
