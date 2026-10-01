"""Bloch-phase plane-wave source for periodic cells (MEEP imported lazily)."""

from __future__ import annotations

import cmath


def bloch_plane_source(cfg):
    """Broadband line source across the full period with transverse phase ``exp(i ky y)``.

    With ``k_point = (0, ky / 2 pi)`` set on the simulation this launches a
    plane wave of transverse wavenumber ``ky`` (exact at the frequency where
    ``ky = 2 pi f sin(theta)``; at other frequencies the same phase front
    radiates the plane wave of the corresponding angle).
    """
    import meep as mp

    freqs = cfg.freqs()
    fcen, df = 0.5 * (freqs[0] + freqs[-1]), freqs[-1] - freqs[0]
    comp = mp.Ez if cfg.polarization == "TM" else mp.Hz
    kwargs = {}
    if cfg.ky != 0:
        ky = cfg.ky
        kwargs["amp_func"] = lambda p: cmath.exp(1j * ky * p.y)
    return [
        mp.Source(
            mp.GaussianSource(fcen, fwidth=cfg.source_fwidth_factor * df),
            component=comp,
            center=mp.Vector3(cfg.source_x, 0),
            size=mp.Vector3(0, cfg.pitch),
            **kwargs,
        )
    ]
