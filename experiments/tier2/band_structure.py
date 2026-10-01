"""Tier 2: angle/Bloch-wavenumber dependence. R/T maps and Harminv ringdown modes vs ky.

For a chosen pitch, scans the transverse wavenumber ky over the first Brillouin zone
(0 .. pi/Lambda) for both polarisations; records R(f, ky), T(f, ky) and the resonant modes
(f, Q) of the structure at each ky.  Results: ``results/tier2/bands``.

    python experiments/tier2/band_structure.py [--pitch 3.0]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from em_rf_lab.modes import find_modes
from em_rf_lab.periodic import PeriodicConfig, run_periodic_or_load

REPO = Path(__file__).resolve().parents[2]


def ky_values(pitch: float, n: int = 11):
    kmax = math.pi / pitch
    return [kmax * i / (n - 1) for i in range(n)]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pitch", type=float, default=3.0)
    ap.add_argument("--resolution", type=float, default=40.0)
    ap.add_argument("--n-ky", type=int, default=11)
    ap.add_argument("--pols", nargs="+", default=["TM"], help="polarisations to scan (TE was not run for the committed results)")
    ap.add_argument("--out", type=Path, default=REPO / "results" / "tier2" / "bands")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    modes_file = args.out / f"modes_L{args.pitch:g}_res{args.resolution:g}.json"
    modes = json.loads(modes_file.read_text()) if modes_file.exists() else {}
    for pol in args.pols:
        for ky in ky_values(args.pitch, args.n_ky):
            cfg = PeriodicConfig(
                name=f"band_{pol}_L{args.pitch:g}_ky{ky:.4f}", polarization=pol, pitch=args.pitch, ky=ky,
                resolution=args.resolution, f_min=0.05, f_max=0.60, n_freq=551, decay_tol=1e-7, max_run_time=1500.0,
            )
            r = run_periodic_or_load(cfg, str(args.out))
            key = f"{pol}_{ky:.4f}"
            if key not in modes:
                ms = find_modes(cfg, fcen=0.30, df=0.50, run_time=1500.0, q_min=3.0)
                modes[key] = [dict(freq=m.freq, q=m.q, amplitude=m.amplitude, error=m.error) for m in ms]
                modes_file.write_text(json.dumps(modes, indent=1))
            print(f"{key}: R/T wall {r.meta['wall_time_total_s']:.0f}s, {len(modes[key])} modes", flush=True)
