"""Matplotlib helpers: consistent, labelled, reproducible scientific figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

__all__ = ["apply_style", "save_fig", "PALETTE"]

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b", "#7f7f7f"]


def apply_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 200,
            "font.size": 10,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.prop_cycle": matplotlib.cycler(color=PALETTE),
            "figure.constrained_layout.use": True,
        }
    )


def save_fig(fig, path: str | Path, caption: str | None = None) -> Path:
    """Save PNG (and a sibling .txt caption file, used by the docs)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    if caption:
        path.with_suffix(".txt").write_text(caption.strip() + "\n")
    plt.close(fig)
    return path
