"""FDTD scattering experiment driver (MEEP).

Workflow for one configuration
------------------------------
1. **Reference run** (no scatterer): records the incident amplitude
   ``F_i(omega)`` at the origin, the incident near-field/flux data on the
   closed monitor box, and the transverse incident profile.
2. **Scattering run** (cylinder present): the stored incident near-field and
   flux data are *subtracted* (``load_minus_*``), so the monitors see only the
   scattered field.  The near-to-far transform of the scattered near field
   gives the far-field echo width; the scattered Poynting flux through the box
   gives the total scattering width independently.
"""

from __future__ import annotations

import math
import time
from typing import Any

import numpy as np

from .config import SimulationConfig
from .farfield import echo_width_from_farfield, extract_axial_farfield, total_width_from_echo
from .geometry import angle_grid, freq_grid, ka_grid
from .io import RunResult, provenance, save_run
from .monitors import attach_monitors
from .sources import plane_wave_sources

__all__ = ["run_scattering", "run_or_load", "meep_available"]


def meep_available() -> bool:
    try:
        import meep  # noqa: F401

        return True
    except Exception:
        return False


class _DecayStop:
    """Stop once the axial field at a probe point has decayed, within run-time bounds.

    Times are measured from the first call (i.e. after the sources switch off,
    because it is used with ``until_after_sources``).
    """

    def __init__(self, cfg: SimulationConfig, component) -> None:
        import meep as mp

        # Decay probe inside the cylinder at radius r*R (default 0.5: the original probe (0.4, 0.3) R).
        # Whispering-gallery modes live near the rim, so r closer to 1 is a stricter decay test.
        r = float(cfg.extra.get("decay_probe_radius", cfg.decay_probe_radius))
        theta = math.atan2(0.3, 0.4)
        probe = mp.Vector3(r * cfg.radius * math.cos(theta), r * cfg.radius * math.sin(theta))
        self._decay = mp.stop_when_fields_decayed(dt=10.0, c=component, pt=probe, decay_by=cfg.decay_tol)
        self._cfg = cfg
        self._t0: float | None = None
        self.hit_max_time = False

    def __call__(self, sim) -> bool:
        decayed = self._decay(sim)  # always call: the closure tracks its own state
        t = sim.meep_time()
        if self._t0 is None:
            self._t0 = t
        elapsed = t - self._t0
        if elapsed < self._cfg.min_run_time:
            return False
        if elapsed >= self._cfg.max_run_time:
            self.hit_max_time = True
            return True
        return bool(decayed)


def _make_sim(cfg: SimulationConfig, with_cylinder: bool):
    import meep as mp

    geometry = []
    if with_cylinder:
        geometry.append(
            mp.Cylinder(radius=cfg.radius, height=mp.inf, center=mp.Vector3(), material=mp.Medium(epsilon=cfg.eps_r))
        )
    symmetries = []
    if cfg.use_symmetry:
        # Ez is even under y -> -y; Hz (a pseudovector component) is odd.
        symmetries.append(mp.Mirror(mp.Y, phase=+1 if cfg.polarization == "TM" else -1))
    return mp.Simulation(
        cell_size=mp.Vector3(cfg.cell_x, cfg.cell_y, 0),
        boundary_layers=[mp.PML(cfg.pml)],
        geometry=geometry,
        sources=plane_wave_sources(cfg),
        resolution=cfg.resolution,
        Courant=cfg.courant,
        symmetries=symmetries,
    )


def _scalar(arr: Any) -> complex:
    return complex(np.asarray(arr).ravel()[0])


def run_scattering(
    cfg: SimulationConfig,
    map_ka: float | None = None,
    save_root: str | None = None,
    verbose: bool = False,
) -> RunResult:
    """Run the reference + scattering simulations and post-process.

    Parameters
    ----------
    cfg : experiment configuration.
    map_ka : if given, also store the complex axial field on a square patch at
        the grid frequency nearest to this size parameter (for near-field
        validation).
    save_root : if given, the result directory is written below this path.

    Returns
    -------
    RunResult with arrays ``ka, freq, angles, F_inc, echo_width, sigma_s_n2f,
    sigma_s_flux, incident_profile_y, profile_y_coords, flux_norm_check`` (and
    ``map_*`` when requested) and metadata including wall times.
    """
    import meep as mp

    mp.verbosity(1 if verbose else 0)
    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    angles = angle_grid(cfg.n_angles)
    ka = ka_grid(cfg)
    freq = freq_grid(cfg)

    map_freq = None
    map_index = None
    if map_ka is not None:
        map_index = int(np.argmin(np.abs(ka - map_ka)))
        map_freq = float(freq[map_index])

    # ---------------------------------------------------------- reference run
    t_wall = time.perf_counter()
    sim_e = _make_sim(cfg, with_cylinder=False)
    mon_e = attach_monitors(sim_e, cfg, map_freq=None)
    # Extra incident-flux check through a transverse line at x = 0.
    import meep as mp  # noqa: F811  (explicit for clarity)
    from .sources import source_band

    fcen, df, _ = source_band(cfg)
    flux_line = sim_e.add_flux(
        fcen, df, cfg.n_freq,
        mp.FluxRegion(center=mp.Vector3(0, 0), size=mp.Vector3(0, 2 * cfg.n2f_half_width), weight=+1),
    )
    stop_e = _DecayStop(cfg, comp)
    sim_e.run(until_after_sources=stop_e)
    t_ref = time.perf_counter() - t_wall
    t_sim_ref = sim_e.meep_time()

    f_inc = np.array([_scalar(sim_e.get_dft_array(mon_e.probe_origin, comp, i)) for i in range(cfg.n_freq)])
    prof_inc = np.array([np.asarray(sim_e.get_dft_array(mon_e.profile_y, comp, i)) for i in range(cfg.n_freq)])
    _, prof_y, _, _ = sim_e.get_array_metadata(dft_cell=mon_e.profile_y)
    n2f_data = sim_e.get_near2far_data(mon_e.n2f)
    flux_data = sim_e.get_flux_data(mon_e.flux)
    flux_line_vals = np.array(mp.get_fluxes(flux_line))
    i0 = 0.5 * np.abs(f_inc) ** 2
    # Calibrate MEEP's flux normalisation: the power through the transverse line must equal
    # the integral of the local plane-wave intensity 0.5|F(y)|^2 measured on the same line.
    y_prof = np.asarray(prof_y)
    sel = np.abs(y_prof) <= cfg.n2f_half_width
    expected = np.trapezoid(0.5 * np.abs(prof_inc[:, sel]) ** 2, y_prof[sel], axis=1)
    # MEEP's DFT flux is Re(E x H*) (no factor 1/2), i.e. twice the time-averaged power, so the
    # calibrated ratio below should be ~1 (it deviates only by discretisation error).
    flux_norm_check = flux_line_vals / (2.0 * expected)

    # ------------------------------------------------------ scattering run
    t_wall = time.perf_counter()
    sim_s = _make_sim(cfg, with_cylinder=True)
    mon_s = attach_monitors(sim_s, cfg, map_freq=map_freq)
    sim_s.load_minus_near2far_data(mon_s.n2f, n2f_data)
    sim_s.load_minus_flux_data(mon_s.flux, flux_data)
    stop_s = _DecayStop(cfg, comp)
    sim_s.run(until_after_sources=stop_s)
    t_scat = time.perf_counter() - t_wall
    t_sim_scat = sim_s.meep_time()

    # ------------------------------------------------------- post-processing
    f_scat = extract_axial_farfield(sim_s, mon_s.n2f, angles, cfg.farfield_distance, cfg.polarization, cfg.n_freq)
    w = echo_width_from_farfield(f_scat, f_inc, cfg.farfield_distance)
    sigma_n2f = total_width_from_echo(w, angles)
    p_scat = np.array(mp.get_fluxes(mon_s.flux))
    sigma_flux = p_scat / (2.0 * i0)  # factor 2: MEEP flux convention (see flux_norm_check)

    arrays: dict[str, np.ndarray] = {
        "ka": ka,
        "freq": freq,
        "angles": angles,
        "F_inc": f_inc,
        "F_scat_far": f_scat,
        "echo_width": w,
        "sigma_s_n2f": sigma_n2f,
        "sigma_s_flux": sigma_flux,
        "incident_profile_y": prof_inc,
        "profile_y_coords": np.asarray(prof_y),
        "flux_norm_check": flux_norm_check,
    }
    if map_freq is not None and mon_s.field_map is not None:
        mx, my, _, _ = sim_s.get_array_metadata(dft_cell=mon_s.field_map)
        fmap = np.asarray(sim_s.get_dft_array(mon_s.field_map, comp, 0))
        arrays["map_x"] = np.asarray(mx)
        arrays["map_y"] = np.asarray(my)
        arrays["map_F"] = fmap / f_inc[map_index]
        arrays["map_ka"] = np.array([ka[map_index]])

    dt = cfg.courant / cfg.resolution
    meta = {
        "wall_time_reference_s": t_ref,
        "wall_time_scattering_s": t_scat,
        "wall_time_total_s": t_ref + t_scat,
        "sim_time_reference": t_sim_ref,
        "sim_time_scattering": t_sim_scat,
        "n_steps_scattering": int(round(t_sim_scat / dt)),
        "dt": dt,
        "decay_reached_reference": not stop_e.hit_max_time,
        "decay_reached_scattering": not stop_s.hit_max_time,
        "cells": int(round(cfg.cell_x * cfg.resolution) * round(cfg.cell_y * cfg.resolution)),
    }
    result = RunResult(config=cfg, arrays=arrays, meta={**provenance(), **meta})
    if save_root is not None:
        save_run(save_root, cfg, arrays, meta)
    return result


def run_or_load(
    cfg: SimulationConfig,
    root: str,
    map_ka: float | None = None,
    force: bool = False,
    verbose: bool = False,
) -> RunResult:
    """Return the cached result for ``cfg`` under ``root`` or run and cache it.

    The cache key is the configuration hash (see :meth:`SimulationConfig.config_hash`),
    so any change to a physical/numerical parameter triggers a fresh run.
    """
    from pathlib import Path

    from .io import load_run, run_dir_name

    path = Path(root) / run_dir_name(cfg)
    if not force and (path / "arrays.npz").exists():
        cached = load_run(path)
        if map_ka is None or "map_F" in cached.arrays:
            return cached
    return run_scattering(cfg, map_ka=map_ka, save_root=root, verbose=verbose)
