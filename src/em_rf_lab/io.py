"""Saving / loading experiment results with full provenance metadata.

A *run directory* contains::

    config.json     the :class:`SimulationConfig` that produced it
    meta.json       provenance: timestamps, code version, library versions, wall time
    arrays.npz      numerical outputs (frequency grid, angles, far fields, ...)

Raw MEEP HDF5 dumps are intentionally *not* stored; only the compact arrays
needed to regenerate every figure are kept.
"""

from __future__ import annotations

import datetime as _dt
import json
import platform
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .config import SimulationConfig, save_config

__all__ = ["RunResult", "provenance", "save_run", "load_run", "run_dir_name"]


@dataclass
class RunResult:
    """Container for one completed experiment."""

    config: Any
    arrays: dict[str, np.ndarray]
    meta: dict[str, Any]

    def __getitem__(self, key: str) -> np.ndarray:
        return self.arrays[key]


def _git_commit(repo: Path | None = None) -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repo or Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return out.stdout.strip() or None if out.returncode == 0 else None
    except Exception:  # pragma: no cover - git absent
        return None


def _git_dirty(repo: Path | None = None) -> bool | None:
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo or Path(__file__).resolve().parent,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return bool(out.stdout.strip()) if out.returncode == 0 else None
    except Exception:  # pragma: no cover
        return None


def provenance(extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Collect environment/version information for reproducibility."""
    import scipy

    info: dict[str, Any] = {
        "timestamp_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "git_commit": _git_commit(),
        "git_dirty": _git_dirty(),
    }
    try:
        import meep

        info["meep"] = getattr(meep, "__version__", "unknown")
    except Exception:
        info["meep"] = None
    if extra:
        info.update(extra)
    return info


def run_dir_name(cfg: Any) -> str:
    """Deterministic directory name: ``<name>_<hash>``."""
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in cfg.name)
    return f"{safe}_{cfg.config_hash()}"


def save_run(root: str | Path, cfg: Any, arrays: dict[str, np.ndarray], meta: dict[str, Any] | None = None) -> Path:
    """Write a run directory under ``root`` and return its path."""
    path = Path(root) / run_dir_name(cfg)
    path.mkdir(parents=True, exist_ok=True)
    save_config(cfg, path / "config.json")
    np.savez_compressed(path / "arrays.npz", **{k: np.asarray(v) for k, v in arrays.items()})
    full_meta = provenance(meta)
    (path / "meta.json").write_text(json.dumps(full_meta, indent=2, sort_keys=True, default=str) + "\n")
    return path


def load_run(path: str | Path, config_cls: Any = SimulationConfig) -> RunResult:
    """Load a run directory written by :func:`save_run`.

    ``config_cls`` must provide ``from_dict`` (default: :class:`SimulationConfig`).
    """
    path = Path(path)
    for name in ("config.json", "meta.json", "arrays.npz"):
        if not (path / name).exists():
            raise FileNotFoundError(f"{path} is missing {name}")
    cfg = config_cls.from_dict(json.loads((path / "config.json").read_text()))
    meta = json.loads((path / "meta.json").read_text())
    with np.load(path / "arrays.npz") as z:
        arrays = {k: z[k] for k in z.files}
    return RunResult(config=cfg, arrays=arrays, meta=meta)
