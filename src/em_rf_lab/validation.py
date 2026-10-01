"""Compare an FDTD run against the exact cylinder (Mie) reference."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import analytical as A
from . import metrics as M
from .io import RunResult

__all__ = ["MieComparison", "compare_to_mie", "mie_echo_width_grid", "mie_sigma_grid"]


def mie_echo_width_grid(ka, angles, radius=1.0, eps_r=4.0, pol="TM") -> np.ndarray:
    """Exact echo width on a (ka, angle) grid, shape ``(n_ka, n_angles)``."""
    return np.array([A.echo_width(angles, k, radius, eps_r, pol) for k in ka])


def mie_sigma_grid(ka, radius=1.0, eps_r=4.0, pol="TM") -> np.ndarray:
    return np.array([A.total_scattering_width(k, radius, eps_r, pol) for k in ka])


@dataclass
class MieComparison:
    ka: np.ndarray
    angles: np.ndarray
    w_ref: np.ndarray               # exact echo width (n_ka, n_ang)
    sigma_ref: np.ndarray           # exact total scattering width (n_ka,)
    rel_l2_angle: np.ndarray        # per-ka relative L2 error of W(phi)
    max_rel_angle: np.ndarray       # per-ka peak-normalised max deviation of W(phi)
    sigma_err_n2f: np.ndarray       # signed relative error of sigma_s from near-to-far
    sigma_err_flux: np.ndarray      # signed relative error of sigma_s from scattered flux
    rel_l2_spectrum: float          # relative L2 error of sigma_s(ka) over the band (n2f)
    rms_rel_l2_angle: float         # RMS over ka of rel_l2_angle
    max_rel_l2_angle: float

    def summary(self) -> dict:
        return {
            "rms_rel_l2_angle": self.rms_rel_l2_angle,
            "max_rel_l2_angle": self.max_rel_l2_angle,
            "rel_l2_spectrum_n2f": self.rel_l2_spectrum,
            "max_abs_sigma_err_n2f": float(np.max(np.abs(self.sigma_err_n2f))),
            "max_abs_sigma_err_flux": float(np.max(np.abs(self.sigma_err_flux))),
        }


def compare_to_mie(run: RunResult) -> MieComparison:
    cfg = run.config
    ka, ang = run["ka"], run["angles"]
    w_ref = mie_echo_width_grid(ka, ang, cfg.radius, cfg.eps_r, cfg.polarization)
    sig_ref = mie_sigma_grid(ka, cfg.radius, cfg.eps_r, cfg.polarization)
    w = run["echo_width"]
    l2 = np.array([M.rel_l2_error(w[i], w_ref[i], ang) for i in range(len(ka))])
    mx = np.array([M.max_rel_deviation(w[i], w_ref[i]) for i in range(len(ka))])
    return MieComparison(
        ka=ka, angles=ang, w_ref=w_ref, sigma_ref=sig_ref,
        rel_l2_angle=l2, max_rel_angle=mx,
        sigma_err_n2f=run["sigma_s_n2f"] / sig_ref - 1.0,
        sigma_err_flux=run["sigma_s_flux"] / sig_ref - 1.0,
        rel_l2_spectrum=M.rel_l2_error(run["sigma_s_n2f"], sig_ref, ka),
        rms_rel_l2_angle=float(np.sqrt(np.mean(l2**2))),
        max_rel_l2_angle=float(np.max(l2)),
    )


# ---------------------------------------------------------------------------
# Numerical-dispersion analysis: rigid frequency shift of the FDTD spectrum
# ---------------------------------------------------------------------------
@dataclass
class ShiftFit:
    delta: float                 # FDTD(ka) ~ Mie(ka * (1 + delta))
    rms_l2_uncompensated: float  # RMS over ka of the relative L2 error of W(phi), as in MieComparison
    rms_l2_compensated: float    # same, after evaluating the exact solution at ka*(1+delta)
    spectrum_l2_uncompensated: float
    spectrum_l2_compensated: float


def fit_frequency_shift(run: RunResult, bound: float = 0.03) -> ShiftFit:
    """Fit a single fractional frequency shift ``delta`` between FDTD and exact spectra.

    Numerical dispersion of the Yee scheme slows the numerical phase velocity, which displaces
    resonances by a (nearly) rigid fraction of frequency, ``delta ~ O(h^2)``.  Because a narrow
    resonance turns a tiny shift into a large local *amplitude* error, separating the shift from the
    residual shape error explains pre-asymptotic convergence near high-Q features.  ``delta`` is
    chosen to minimise the spectrum error of the total scattering width (near-to-far estimate).
    """
    from scipy.optimize import minimize_scalar

    cfg = run.config
    ka, ang = run["ka"], run["angles"]
    sig = run["sigma_s_n2f"]

    def spec_err(d: float) -> float:
        ref = mie_sigma_grid(ka * (1 + d), cfg.radius, cfg.eps_r, cfg.polarization)
        return M.rel_l2_error(sig, ref, ka)

    sol = minimize_scalar(spec_err, bounds=(-bound, bound), method="bounded", options={"xatol": 1e-7})
    d = float(sol.x)

    def rms_angle(dd: float) -> float:
        ref = mie_echo_width_grid(ka * (1 + dd), ang, cfg.radius, cfg.eps_r, cfg.polarization)
        e = [M.rel_l2_error(run["echo_width"][i], ref[i], ang) for i in range(len(ka))]
        return float(np.sqrt(np.mean(np.square(e))))

    return ShiftFit(d, rms_angle(0.0), rms_angle(d), spec_err(0.0), float(sol.fun))
