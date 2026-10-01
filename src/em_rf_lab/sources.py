"""Source construction (MEEP imported lazily so the package imports without it)."""

from __future__ import annotations

from .config import SimulationConfig
from .geometry import freq_grid

__all__ = ["source_band", "plane_wave_sources"]


def source_band(cfg: SimulationConfig) -> tuple[float, float, float]:
    """Return ``(fcen, df, fwidth)`` for the Gaussian pulse.

    ``fcen``/``df`` describe the monitored frequency band (centre and full
    width); ``fwidth`` is the Gaussian-pulse bandwidth parameter passed to MEEP.
    """
    f = freq_grid(cfg)
    fmin, fmax = float(f[0]), float(f[-1])
    fcen = 0.5 * (fmin + fmax)
    df = fmax - fmin
    fwidth = cfg.source_fwidth_factor * df if df > 0 else 0.2 * fcen
    return fcen, df, fwidth


def plane_wave_sources(cfg: SimulationConfig):
    """Broadband line source radiating a quasi-plane wave toward +x.

    The source spans the full cell height (into the PML).  By default its
    amplitude is Gaussian in ``y`` (1/e half-width ``beam_waist``), which removes
    the aperture-edge diffraction of a hard-edged line source; the residual
    departure from a plane wave over the scatterer is ``O((R/beam_waist)^2)``
    and is quantified in the domain-size study.  ``beam_waist = 0`` selects a
    uniform (hard-edged) line.  TM uses ``Ez``; TE uses ``Hz``.
    """
    import meep as mp

    import math

    fcen, _, fwidth = source_band(cfg)
    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    kwargs = {}
    if cfg.beam_waist > 0:
        w = cfg.beam_waist
        kwargs["amp_func"] = lambda p: math.exp(-((p.y / w) ** 2))
    return [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=fwidth),
            component=comp,
            center=mp.Vector3(cfg.source_x, 0, 0),
            size=mp.Vector3(0, cfg.cell_y, 0),
            **kwargs,
        )
    ]
