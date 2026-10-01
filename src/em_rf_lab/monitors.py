"""Monitor construction: near-to-far surface, scattered-power flux box, field probes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .config import SimulationConfig
from .geometry import near2far_segments
from .sources import source_band

__all__ = ["MonitorSet", "attach_monitors"]


@dataclass
class MonitorSet:
    """Handles to the monitors attached to one simulation."""

    n2f: Any
    flux: Any
    probe_origin: Any          # DFT of the axial field at the origin (incident amplitude)
    profile_y: Any             # DFT of the axial field along the line x = 0 (plane-wave quality)
    field_map: Any | None      # DFT of the axial field on a square patch at one frequency
    field_map_freq: float | None


def attach_monitors(sim, cfg: SimulationConfig, map_freq: float | None = None) -> MonitorSet:
    """Attach all DFT monitors to ``sim`` (a :class:`meep.Simulation`)."""
    import meep as mp

    fcen, df, _ = source_band(cfg)
    nfreq = cfg.n_freq
    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    segs = near2far_segments(cfg.n2f_half_width)

    n2f = sim.add_near2far(
        fcen, df, nfreq,
        *[mp.Near2FarRegion(center=mp.Vector3(*s["center"]), size=mp.Vector3(*s["size"]), weight=s["weight"]) for s in segs],
    )
    flux = sim.add_flux(
        fcen, df, nfreq,
        *[mp.FluxRegion(center=mp.Vector3(*s["center"]), size=mp.Vector3(*s["size"]), weight=s["weight"]) for s in segs],
    )
    probe_origin = sim.add_dft_fields([comp], fcen, df, nfreq, where=mp.Volume(center=mp.Vector3(), size=mp.Vector3()))

    inner_y = cfg.cell_y / 2 - cfg.pml
    profile_y = sim.add_dft_fields(
        [comp], fcen, df, nfreq,
        where=mp.Volume(center=mp.Vector3(0, 0), size=mp.Vector3(0, 2 * inner_y)),
    )

    field_map = None
    if map_freq is not None:
        half = cfg.n2f_half_width
        field_map = sim.add_dft_fields(
            [comp], map_freq, 0, 1,
            where=mp.Volume(center=mp.Vector3(), size=mp.Vector3(2 * half, 2 * half)),
        )
    return MonitorSet(n2f, flux, probe_origin, profile_y, field_map, map_freq)
