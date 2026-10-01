"""Validated, hashable configuration for a scattering experiment.

All lengths are in MEEP length units (the cylinder radius sets the scale; with
``radius = 1`` lengths are in units of the cylinder radius).  Frequencies are
in MEEP units ``c/L``; with ``c = 1`` the vacuum wavenumber is ``k = 2 pi f``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["SimulationConfig", "ConfigError", "load_config", "save_config"]


class ConfigError(ValueError):
    """Raised for physically or numerically invalid configurations."""


@dataclass(frozen=True)
class SimulationConfig:
    """Parameters of one FDTD scattering experiment (dielectric cylinder).

    Attributes
    ----------
    name : human-readable label (not part of the physics hash).
    radius, eps_r : cylinder radius and (real, lossless) relative permittivity.
    polarization : ``"TM"`` (E_z source) or ``"TE"`` (H_z source).
    resolution : FDTD pixels per MEEP length unit.
    ka_min, ka_max, n_freq : uniformly spaced size-parameter grid ``k R``.
    cell_x, cell_y : full cell size *including* PML.
    pml : PML thickness (all four sides).
    source_gap : distance between the inner PML edge and the source line.
    beam_waist : 1/e amplitude half-width of the Gaussian source profile along y
        (0 selects a uniform, hard-edged line source).
    n2f_half_width : half side length of the near-to-far box (centred on origin).
    farfield_distance : radius at which far fields are evaluated (units of length).
    decay_tol : field-decay tolerance used as the run-termination criterion.
    max_run_time, min_run_time : bounds (MEEP time units, measured after the source
        switches off) around the decay test.  Narrow high-Q resonances ring for
        ~Q/(pi f) time units, so ``min_run_time`` must exceed several such times
        (see the verification study of run-time convergence).
    decay_probe_radius : radius (units of the cylinder radius) of the interior
        point used by the decay test; whispering-gallery modes live near the rim,
        so values near 1 are the stricter test.
    source_fwidth_factor : Gaussian source bandwidth as a multiple of the band width.
    n_angles : number of observation angles in ``[0, 2 pi)``.
    courant : Courant factor (MEEP default 0.5).
    use_symmetry : exploit the y -> -y mirror symmetry of the problem (halves the cost).
    """

    name: str = "baseline"
    radius: float = 1.0
    eps_r: float = 4.0
    polarization: str = "TM"
    resolution: float = 30.0
    ka_min: float = 2.0
    ka_max: float = 6.0
    n_freq: int = 41
    cell_x: float = 20.0
    cell_y: float = 56.0
    pml: float = 5.0
    source_gap: float = 0.5
    beam_waist: float = 20.0
    n2f_half_width: float = 2.5
    farfield_distance: float = 1.0e5
    decay_tol: float = 1.0e-6
    max_run_time: float = 3000.0
    min_run_time: float = 500.0
    decay_probe_radius: float = 0.9
    source_fwidth_factor: float = 1.5
    n_angles: int = 91
    courant: float = 0.5
    use_symmetry: bool = True
    extra: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ checks
    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        def need(cond: bool, msg: str) -> None:
            if not cond:
                raise ConfigError(msg)

        need(self.polarization in ("TM", "TE"), f"polarization must be TM or TE, got {self.polarization!r}")
        need(self.radius > 0, "radius must be positive")
        need(self.eps_r >= 1.0, "eps_r must be >= 1 (lossless dielectric)")
        need(self.resolution >= 4, "resolution must be >= 4 pixels per unit length")
        need(0 < self.ka_min < self.ka_max, "require 0 < ka_min < ka_max")
        need(self.n_freq >= 1, "n_freq must be >= 1")
        need(self.pml > 0, "pml thickness must be positive")
        need(self.cell_x > 2 * self.pml and self.cell_y > 2 * self.pml, "cell must be larger than twice the PML thickness")
        need(self.source_gap > 0, "source_gap must be positive")
        need(self.beam_waist >= 0, "beam_waist must be >= 0 (0 = uniform line source)")
        need(self.n2f_half_width > self.radius, "near-to-far box must enclose the cylinder")
        inner_x = self.cell_x / 2 - self.pml
        inner_y = self.cell_y / 2 - self.pml
        need(self.n2f_half_width < inner_y, "near-to-far box must lie inside the PML-free region (y)")
        src_x = self.source_x
        need(-inner_x < src_x < -self.n2f_half_width, "source line must lie between the PML and the near-to-far box")
        need(self.farfield_distance > 100 * self.radius, "farfield_distance should be >> radius (>=100 R)")
        need(self.decay_tol > 0, "decay_tol must be positive")
        need(0 < self.decay_probe_radius < 1.0, "decay_probe_radius must lie in (0, 1)")
        need(0 < self.courant <= 0.5, "courant must be in (0, 0.5]")
        need(self.max_run_time > self.min_run_time > 0, "require 0 < min_run_time < max_run_time")
        need(self.n_angles >= 3, "n_angles must be >= 3")
        need(self.source_fwidth_factor > 0, "source_fwidth_factor must be positive")

    # ---------------------------------------------------------------- derived
    @property
    def source_x(self) -> float:
        """x position of the planar (line) source."""
        return -(self.cell_x / 2 - self.pml - self.source_gap)

    @property
    def components(self) -> str:
        """Axial field component used as source / observable."""
        return "Ez" if self.polarization == "TM" else "Hz"

    # ------------------------------------------------------------ serialisation
    def to_dict(self) -> dict[str, Any]:
        return dataclasses.asdict(self)

    def physics_dict(self) -> dict[str, Any]:
        """Dictionary used for hashing (excludes the free-text ``name``)."""
        d = self.to_dict()
        d.pop("name", None)
        return d

    def config_hash(self, length: int = 10) -> str:
        blob = json.dumps(self.physics_dict(), sort_keys=True, default=float).encode()
        return hashlib.sha256(blob).hexdigest()[:length]

    def replace(self, **changes: Any) -> "SimulationConfig":
        return dataclasses.replace(self, **changes)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "SimulationConfig":
        names = {f.name for f in dataclasses.fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ConfigError(f"unknown configuration keys: {sorted(unknown)}")
        return cls(**d)


def save_config(cfg: SimulationConfig, path: str | Path) -> None:
    Path(path).write_text(json.dumps(cfg.to_dict(), indent=2, sort_keys=True) + "\n")


def load_config(path: str | Path) -> SimulationConfig:
    return SimulationConfig.from_dict(json.loads(Path(path).read_text()))
