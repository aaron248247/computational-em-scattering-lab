"""Resonant modes of a periodic layer from time-domain ringdown (Harminv).

A broadband point source excites the structure at a fixed Bloch wavenumber
``ky``; after the source has decayed, Harminv fits the field ringdown at probe
points with a sum of damped exponentials, returning complex frequencies
``f_r, Q = f_r / (2 |Im f|)``.  This is independent of the reflectance
line-shape fits in :mod:`em_rf_lab.resonance`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .periodic import PeriodicConfig

__all__ = ["Mode", "find_modes"]


@dataclass(frozen=True)
class Mode:
    freq: float
    q: float
    amplitude: float
    error: float


def find_modes(
    cfg: PeriodicConfig,
    fcen: float,
    df: float,
    run_time: float = 600.0,
    q_min: float = 5.0,
    err_max: float = 1e-3,
) -> list[Mode]:
    """Modes of the *structure* (no incident plane wave) near ``fcen`` within ``df``.

    Probe/source points are placed at generic (non-symmetric) positions inside
    the unit cell so that modes of every symmetry class are excited.
    """
    import meep as mp

    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    eps = mp.Medium(epsilon=cfg.eps_r)
    if cfg.structure == "cylinder":
        geometry = [mp.Cylinder(radius=cfg.radius, height=mp.inf, material=eps)]
    else:
        geometry = [mp.Block(size=mp.Vector3(cfg.thickness, mp.inf, mp.inf), material=eps)]
    src_pts = [mp.Vector3(0.17 * cfg.half_extent, 0.31 * cfg.pitch), mp.Vector3(-0.23 * cfg.half_extent, -0.19 * cfg.pitch)]
    probe_pts = [mp.Vector3(0.29 * cfg.half_extent, -0.13 * cfg.pitch), mp.Vector3(-0.11 * cfg.half_extent, 0.37 * cfg.pitch)]
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=comp, center=p) for p in src_pts]
    kp = mp.Vector3(0, cfg.ky / (2 * math.pi), 0) if cfg.ky != 0 else mp.Vector3()
    sim = mp.Simulation(
        cell_size=mp.Vector3(cfg.cell_x, cfg.pitch, 0),
        boundary_layers=[mp.PML(cfg.pml, direction=mp.X)],
        geometry=geometry,
        sources=sources,
        k_point=kp,
        resolution=cfg.resolution,
        Courant=cfg.courant,
    )
    mp.verbosity(0)
    hs = [mp.Harminv(comp, p, fcen, df) for p in probe_pts]
    sim.run(mp.after_sources(*hs), until_after_sources=run_time)
    found: dict[int, Mode] = {}
    for h in hs:
        for m in h.modes:
            if m.Q < q_min or abs(m.err) > err_max * abs(m.freq) or abs(m.freq - fcen) > df / 2:
                continue
            # harminv also returns negative-Q (growing) fits: drop them
            mode = Mode(float(m.freq), float(m.Q), float(abs(m.amp)), float(abs(m.err)))
            key = int(round(mode.freq / 2e-4))
            if key not in found or mode.amplitude > found[key].amplitude:
                found[key] = mode
    return sorted(found.values(), key=lambda m: m.freq)
