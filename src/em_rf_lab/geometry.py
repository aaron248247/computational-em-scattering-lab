"""Pure-Python geometry and unit conversions (no MEEP import).

Unit system: MEEP units with ``c = eps_0 = mu_0 = 1``.  A length ``L`` is in
whatever unit the user picks (the cylinder radius by default); frequency is
``f = omega / 2 pi`` in units of ``c / L`` and the vacuum wavenumber is
``k = 2 pi f``.
"""

from __future__ import annotations

import math

import numpy as np

from .config import SimulationConfig

__all__ = [
    "ka_to_freq",
    "freq_to_ka",
    "vacuum_wavelength",
    "points_per_wavelength",
    "ka_grid",
    "freq_grid",
    "angle_grid",
    "near2far_segments",
    "cell_size",
    "courant_timestep",
]


def ka_to_freq(ka: float | np.ndarray, radius: float = 1.0):
    """Size parameter ``k R`` -> MEEP frequency ``f = k / (2 pi) = ka / (2 pi R)``."""
    return np.asarray(ka) / (2.0 * math.pi * radius) if np.ndim(ka) else ka / (2.0 * math.pi * radius)


def freq_to_ka(freq: float | np.ndarray, radius: float = 1.0):
    """MEEP frequency -> size parameter ``k R = 2 pi f R``."""
    return 2.0 * math.pi * np.asarray(freq) * radius if np.ndim(freq) else 2.0 * math.pi * freq * radius


def vacuum_wavelength(freq: float) -> float:
    """``lambda = c / f`` (c = 1)."""
    if freq <= 0:
        raise ValueError("frequency must be positive")
    return 1.0 / freq


def points_per_wavelength(resolution: float, freq: float, refractive_index: float = 1.0) -> float:
    """Pixels per *in-medium* wavelength ``lambda / n``."""
    return resolution * vacuum_wavelength(freq) / refractive_index


def ka_grid(cfg: SimulationConfig) -> np.ndarray:
    """Uniform size-parameter grid (inclusive of both ends)."""
    return np.linspace(cfg.ka_min, cfg.ka_max, cfg.n_freq)


def freq_grid(cfg: SimulationConfig) -> np.ndarray:
    return ka_to_freq(ka_grid(cfg), cfg.radius)


def angle_grid(n_angles: int) -> np.ndarray:
    """Observation angles in ``[0, pi]`` (forward to backscatter), ``n`` points.

    The echo width of a circular cylinder illuminated along +x is symmetric
    about the x axis, so ``[0, pi]`` carries all information.
    """
    if n_angles < 3:
        raise ValueError("need at least 3 angles")
    return np.linspace(0.0, math.pi, n_angles)


def cell_size(cfg: SimulationConfig) -> tuple[float, float]:
    return (cfg.cell_x, cfg.cell_y)


def courant_timestep(cfg: SimulationConfig) -> float:
    """FDTD time step ``dt = S / resolution`` (c = 1) in MEEP units."""
    return cfg.courant / cfg.resolution


def near2far_segments(half_width: float) -> list[dict]:
    """The four faces of the square near-to-far box as (center, size, normal-weight).

    Returns dicts with ``center`` and ``size`` (2-tuples), ``direction`` (the
    normal axis, ``"x"`` or ``"y"``) and ``weight`` (+1 for the outward normal
    along +axis, -1 for -axis; MEEP's Near2FarRegion ``weight`` convention).
    """
    h = float(half_width)
    return [
        {"center": (0.0, +h), "size": (2 * h, 0.0), "direction": "y", "weight": +1},
        {"center": (0.0, -h), "size": (2 * h, 0.0), "direction": "y", "weight": -1},
        {"center": (+h, 0.0), "size": (0.0, 2 * h), "direction": "x", "weight": +1},
        {"center": (-h, 0.0), "size": (0.0, 2 * h), "direction": "x", "weight": -1},
    ]
