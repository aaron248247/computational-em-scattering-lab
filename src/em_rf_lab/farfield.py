"""Far-field extraction and echo-width normalisation."""

from __future__ import annotations

import numpy as np

__all__ = ["echo_width_from_farfield", "total_width_from_echo", "extract_axial_farfield"]

# Index of the axial component within MEEP's (Ex, Ey, Ez, Hx, Hy, Hz) far-field vector.
_AXIAL_INDEX = {"TM": 2, "TE": 5}


def echo_width_from_farfield(f_scat: np.ndarray, f_inc: np.ndarray, rho: float) -> np.ndarray:
    """2-D echo width ``W = 2 pi rho |F_s|^2 / |F_i|^2`` (units of ``rho``).

    Parameters
    ----------
    f_scat : complex scattered axial field at distance ``rho``; shape ``(n_freq, n_angles)``.
    f_inc : complex incident axial amplitude at the origin; shape ``(n_freq,)``.
    rho : observation distance (same length unit as the geometry).
    """
    f_scat = np.asarray(f_scat)
    f_inc = np.asarray(f_inc)
    if f_scat.ndim != 2 or f_inc.ndim != 1 or f_scat.shape[0] != f_inc.shape[0]:
        raise ValueError("expected f_scat (n_freq, n_angles) and f_inc (n_freq,)")
    if rho <= 0:
        raise ValueError("rho must be positive")
    return 2.0 * np.pi * rho * np.abs(f_scat) ** 2 / np.abs(f_inc)[:, None] ** 2


def total_width_from_echo(w: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """Total scattering width ``sigma_s = (1/2pi) \\oint W dphi`` from ``W`` on ``[0, pi]``.

    Uses the up-down symmetry of a cylinder illuminated along +x
    (``W(-phi) = W(phi)``): ``sigma_s = (1/pi) \\int_0^pi W dphi`` (trapezoid
    rule, which is spectrally accurate for the even periodic extension).
    """
    w = np.asarray(w, dtype=float)
    angles = np.asarray(angles, dtype=float)
    if angles[0] != 0.0 or not np.isclose(angles[-1], np.pi):
        raise ValueError("angles must span [0, pi]")
    return np.trapezoid(w, angles, axis=-1) / np.pi


def extract_axial_farfield(sim, n2f, angles: np.ndarray, rho: float, polarization: str, n_freq: int) -> np.ndarray:
    """Evaluate the axial scattered far field at ``rho`` on ``angles``.

    Returns a complex array ``(n_freq, n_angles)``.  ``sim`` must already have
    the incident near-field data subtracted (``load_minus_near2far_data``).
    """
    import meep as mp

    idx = _AXIAL_INDEX[polarization]
    out = np.empty((n_freq, len(angles)), dtype=complex)
    for j, phi in enumerate(angles):
        ff = np.asarray(sim.get_farfield(n2f, mp.Vector3(rho * np.cos(phi), rho * np.sin(phi))))
        out[:, j] = ff.reshape(n_freq, 6)[:, idx]
    return out
