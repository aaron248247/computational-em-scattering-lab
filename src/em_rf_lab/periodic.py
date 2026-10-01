"""Tier 2: one-dimensional periodic layers (frequency-selective surface / metasurface).

Geometry (2-D, invariant along z)
---------------------------------
A layer periodic in ``y`` with pitch ``Lambda``; light propagates along ``+x``.
Two structures are supported:

* ``"cylinder"`` - an infinite array of circular dielectric cylinders (radius
  ``R``, permittivity ``eps_r``; *the same scatterer that is verified against
  exact Mie theory in Tier 1*) centred on ``x = 0``.
* ``"slab"`` - a homogeneous dielectric film of thickness ``d``; used only as a
  verification structure because it has an exact transfer-matrix solution.

Polarisation naming follows Tier 1: ``"TM"`` = field ``E_z`` (E parallel to the
cylinder axis, i.e. *s*-polarised for a plane of incidence normal to z) and
``"TE"`` = ``H_z`` (*p*-polarised).

Observables
-----------
Total reflectance ``R(f)`` and transmittance ``T(f)`` (summed over all
propagating diffraction orders) from scattered/transmitted Poynting flux divided
by the incident flux of a reference (empty-cell) run; ``A = 1 - R - T`` (should
vanish for a lossless layer).
"""

from __future__ import annotations

import cmath
import dataclasses
import hashlib
import json
import math
import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import ConfigError
from .io import RunResult, provenance, save_run

__all__ = [
    "PeriodicConfig",
    "slab_rt",
    "diffraction_orders",
    "wood_anomaly_frequencies",
    "incidence_angle",
    "run_periodic",
    "run_periodic_or_load",
]


@dataclass(frozen=True)
class PeriodicConfig:
    """Parameters of one periodic-layer reflectance/transmittance experiment.

    Lengths are in MEEP units (cylinder radius ``radius`` sets the scale in the
    benchmark-consistent studies); ``ky`` is the transverse (Bloch) wavenumber
    ``k_y = 2 pi f sin(theta)`` in radians per unit length, so ``ky = 0`` is
    normal incidence.  Because ``ky`` is held fixed while frequency is swept, the
    incidence angle varies with frequency (``sin theta = ky / (2 pi f)``).
    """

    name: str = "array"
    structure: str = "cylinder"           # "cylinder" | "slab"
    pitch: float = 3.0
    radius: float = 1.0
    thickness: float = 2.0                # slab thickness (structure == "slab")
    eps_r: float = 4.0
    polarization: str = "TM"
    ky: float = 0.0
    resolution: float = 40.0
    f_min: float = 0.05
    f_max: float = 0.60
    n_freq: int = 201
    margin: float = 6.0                   # structure edge -> inner PML edge
    monitor_gap: float = 3.0              # structure edge -> R/T monitor planes
    pml: float = 4.0
    decay_tol: float = 1.0e-7
    max_run_time: float = 3000.0
    min_run_time: float = 20.0
    source_fwidth_factor: float = 1.3
    courant: float = 0.5
    extra: dict[str, Any] = dataclasses.field(default_factory=dict)

    def __post_init__(self) -> None:
        def need(c: bool, m: str) -> None:
            if not c:
                raise ConfigError(m)

        need(self.structure in ("cylinder", "slab"), "structure must be 'cylinder' or 'slab'")
        need(self.polarization in ("TM", "TE"), "polarization must be TM or TE")
        need(self.pitch > 0 and self.eps_r >= 1.0, "pitch must be > 0 and eps_r >= 1")
        need(self.resolution >= 8, "resolution must be >= 8")
        need(0 < self.f_min < self.f_max, "require 0 < f_min < f_max")
        need(self.n_freq >= 2, "n_freq must be >= 2")
        need(self.pml > 0 and self.margin > self.monitor_gap > 0, "require margin > monitor_gap > 0 and pml > 0")
        need(self.margin - self.monitor_gap >= 1.0, "source must sit >= 1 unit beyond the monitor planes")
        need(0 < self.courant <= 0.5, "courant must be in (0, 0.5]")
        need(self.max_run_time > self.min_run_time > 0 and self.decay_tol > 0, "invalid run-time settings")
        if self.structure == "cylinder":
            need(self.radius > 0 and 2 * self.radius < self.pitch, "cylinders must not overlap: 2*radius < pitch")
        else:
            need(self.thickness > 0, "slab thickness must be positive")
        # the highest frequency must be propagating for the requested ky at least somewhere in the band
        need(abs(self.ky) >= 0, "ky must be finite")
        need(math.isfinite(self.ky), "ky must be finite")

    # ---------------------------------------------------------------- derived
    @property
    def half_extent(self) -> float:
        return self.radius if self.structure == "cylinder" else self.thickness / 2

    @property
    def cell_x(self) -> float:
        return 2 * (self.half_extent + self.margin + self.pml)

    @property
    def monitor_x(self) -> float:
        return self.half_extent + self.monitor_gap

    @property
    def source_x(self) -> float:
        return -(self.half_extent + self.margin - 0.5)

    def freqs(self) -> np.ndarray:
        return np.linspace(self.f_min, self.f_max, self.n_freq)

    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)

    def config_hash(self, length: int = 10) -> str:
        d = self.to_dict()
        d.pop("name", None)
        return hashlib.sha256(json.dumps(d, sort_keys=True, default=float).encode()).hexdigest()[:length]

    def replace(self, **kw: Any) -> "PeriodicConfig":
        return dataclasses.replace(self, **kw)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "PeriodicConfig":
        names = {f.name for f in dataclasses.fields(cls)}
        bad = set(d) - names
        if bad:
            raise ConfigError(f"unknown configuration keys: {sorted(bad)}")
        return cls(**d)


# ---------------------------------------------------------------------------
# Exact reference: homogeneous slab (Airy / transfer-matrix), vacuum on both sides
# ---------------------------------------------------------------------------
def incidence_angle(f: float | np.ndarray, ky: float) -> np.ndarray:
    """``theta = asin(ky / (2 pi f))`` (NaN where the wave is evanescent)."""
    s = ky / (2 * np.pi * np.asarray(f, dtype=float))
    with np.errstate(invalid="ignore"):
        return np.where(np.abs(s) <= 1.0, np.arcsin(np.clip(s, -1, 1)), np.nan)


def slab_rt(f, n: float, d: float, ky: float = 0.0, polarization: str = "TM"):
    """Exact ``(R, T)`` of a lossless slab of index ``n`` and thickness ``d`` in vacuum.

    ``polarization="TM"`` (E_z, *s*) or ``"TE"`` (H_z, *p*).  ``ky`` is the conserved
    transverse wavenumber.  Frequencies where the incident wave is evanescent return NaN.
    """
    if polarization not in ("TM", "TE"):
        raise ValueError("polarization must be TM or TE")
    f = np.atleast_1d(np.asarray(f, dtype=float))
    R = np.full(f.shape, np.nan)
    T = np.full(f.shape, np.nan)
    for i, fi in enumerate(f):
        k0 = 2 * np.pi * fi
        if abs(ky) >= k0:
            continue
        kx1 = math.sqrt(k0**2 - ky**2)
        kx2 = cmath.sqrt((n * k0) ** 2 - ky**2)
        if polarization == "TM":        # s: r = (kx1 - kx2)/(kx1 + kx2)
            r12 = (kx1 - kx2) / (kx1 + kx2)
            t12 = 2 * kx1 / (kx1 + kx2)
            t21 = 2 * kx2 / (kx1 + kx2)
        else:                           # p (mu = 1): r = (eps2 kx1 - eps1 kx2)/(eps2 kx1 + eps1 kx2)
            e2 = n**2
            r12 = (e2 * kx1 - kx2) / (e2 * kx1 + kx2)
            t12 = 2 * e2 * kx1 / (e2 * kx1 + kx2) * 1.0
            t21 = 2 * kx2 / (e2 * kx1 + kx2) * 1.0
        r23 = -r12
        ph = cmath.exp(1j * kx2 * d)
        den = 1 + r12 * r23 * ph**2
        r = (r12 + r23 * ph**2) / den
        t = t12 * t21 * ph / den
        # power transmittance: |t|^2 for identical outer media (t defined for the tangential-field
        # amplitude in a vacuum-vacuum sandwich: t12*t21 is already amplitude-to-amplitude)
        R[i] = abs(r) ** 2
        T[i] = abs(t) ** 2
    return R, T


# ---------------------------------------------------------------------------
# Grating bookkeeping
# ---------------------------------------------------------------------------
def diffraction_orders(f: float, pitch: float, ky: float = 0.0) -> list[int]:
    """Indices ``m`` of propagating orders: ``|ky + 2 pi m / pitch| < 2 pi f``."""
    k0 = 2 * np.pi * f
    m_max = int(np.ceil((k0 + abs(ky)) * pitch / (2 * np.pi))) + 1
    return [m for m in range(-m_max, m_max + 1) if abs(ky + 2 * np.pi * m / pitch) < k0]


def wood_anomaly_frequencies(pitch: float, ky: float = 0.0, m_max: int = 3) -> list[tuple[int, float]]:
    """Order-emergence (Rayleigh-Wood) frequencies ``f_m = |ky + 2 pi m / pitch| / (2 pi)`` for ``m != 0``."""
    out = [(m, abs(ky + 2 * np.pi * m / pitch) / (2 * np.pi)) for m in range(-m_max, m_max + 1) if m != 0]
    return sorted(out, key=lambda t: t[1])


# ---------------------------------------------------------------------------
# FDTD
# ---------------------------------------------------------------------------
def _build(cfg: PeriodicConfig, with_structure: bool):
    import meep as mp
    from .periodic_sources import bloch_plane_source  # local to keep this module importable without meep

    geometry = []
    if with_structure:
        eps = mp.Medium(epsilon=cfg.eps_r)
        if cfg.structure == "cylinder":
            geometry = [mp.Cylinder(radius=cfg.radius, height=mp.inf, center=mp.Vector3(), material=eps)]
        else:
            geometry = [mp.Block(size=mp.Vector3(cfg.thickness, mp.inf, mp.inf), center=mp.Vector3(), material=eps)]
    kp = mp.Vector3(0, cfg.ky / (2 * math.pi), 0) if cfg.ky != 0 else mp.Vector3()
    return mp.Simulation(
        cell_size=mp.Vector3(cfg.cell_x, cfg.pitch, 0),
        boundary_layers=[mp.PML(cfg.pml, direction=mp.X)],
        geometry=geometry,
        sources=bloch_plane_source(cfg),
        k_point=kp,
        resolution=cfg.resolution,
        Courant=cfg.courant,
    )


class _Stop:
    def __init__(self, cfg: PeriodicConfig, comp) -> None:
        import meep as mp

        self._d = mp.stop_when_fields_decayed(dt=20.0, c=comp, pt=mp.Vector3(0.1, 0.3 * cfg.pitch), decay_by=cfg.decay_tol)
        self._cfg, self._t0, self.hit_max = cfg, None, False

    def __call__(self, sim) -> bool:
        ok = self._d(sim)
        t = sim.meep_time()
        if self._t0 is None:
            self._t0 = t
        el = t - self._t0
        if el < self._cfg.min_run_time:
            return False
        if el >= self._cfg.max_run_time:
            self.hit_max = True
            return True
        return bool(ok)


def run_periodic(cfg: PeriodicConfig, save_root: str | None = None, verbose: bool = False) -> RunResult:
    """Run reference + structure simulations and return R(f), T(f), A(f)."""
    import meep as mp

    mp.verbosity(1 if verbose else 0)
    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    freqs = cfg.freqs()
    fcen, df = 0.5 * (freqs[0] + freqs[-1]), freqs[-1] - freqs[0]

    def monitors(sim):
        fr = sim.add_flux(fcen, df, cfg.n_freq, mp.FluxRegion(center=mp.Vector3(-cfg.monitor_x), size=mp.Vector3(0, cfg.pitch)))
        ft = sim.add_flux(fcen, df, cfg.n_freq, mp.FluxRegion(center=mp.Vector3(+cfg.monitor_x), size=mp.Vector3(0, cfg.pitch)))
        return fr, ft

    t0 = time.perf_counter()
    s0 = _build(cfg, with_structure=False)
    fr0, ft0 = monitors(s0)
    st0 = _Stop(cfg, comp)
    s0.run(until_after_sources=st0)
    inc_r = np.array(mp.get_fluxes(fr0))
    inc_t = np.array(mp.get_fluxes(ft0))
    data_r = s0.get_flux_data(fr0)
    t_ref = time.perf_counter() - t0

    t0 = time.perf_counter()
    s1 = _build(cfg, with_structure=True)
    fr1, ft1 = monitors(s1)
    s1.load_minus_flux_data(fr1, data_r)
    st1 = _Stop(cfg, comp)
    s1.run(until_after_sources=st1)
    refl = -np.array(mp.get_fluxes(fr1)) / inc_r
    trans = np.array(mp.get_fluxes(ft1)) / inc_t
    t_str = time.perf_counter() - t0

    arrays = {
        "freq": freqs,
        "R": refl,
        "T": trans,
        "A": 1.0 - refl - trans,
        "incident_flux": inc_r,
        "theta": incidence_angle(freqs, cfg.ky),
    }
    meta = {
        "wall_time_total_s": t_ref + t_str,
        "sim_time_structure": s1.meep_time(),
        "decay_reached": not st1.hit_max,
        "cells": int(round(cfg.cell_x * cfg.resolution) * round(cfg.pitch * cfg.resolution)),
    }
    res = RunResult(config=cfg, arrays=arrays, meta={**provenance(), **meta})
    if save_root is not None:
        save_run(save_root, cfg, arrays, meta)
    return res


def run_periodic_or_load(cfg: PeriodicConfig, root: str, force: bool = False, verbose: bool = False) -> RunResult:
    from pathlib import Path

    from .io import load_run, run_dir_name

    path = Path(root) / run_dir_name(cfg)
    if not force and (path / "arrays.npz").exists():
        return load_run(path, PeriodicConfig)
    return run_periodic(cfg, save_root=root, verbose=verbose)
