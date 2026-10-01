"""Exact (Mie-series) scattering of a plane wave by an infinite circular cylinder.

This module is the *independent reference* against which the FDTD results are
validated.  It depends only on NumPy/SciPy and contains no simulation code.

Geometry and conventions
------------------------
* Infinite, homogeneous circular cylinder of radius ``R`` along ``z``; the
  surrounding medium is vacuum (eps = mu = 1).
* Plane wave incident along ``+x`` at normal incidence, time dependence
  ``exp(-i w t)`` (the MEEP convention), so outgoing waves are ``H_n^(1)``.
* Two scalar polarisations, where ``F`` is the field component parallel to
  the cylinder axis:

  - ``"TM"`` : ``F = E_z``  (electric field parallel to the axis)
  - ``"TE"`` : ``F = H_z``  (magnetic field parallel to the axis)

Expansion (per unit incident amplitude)::

    F_inc = exp(i k x)            = sum_n i^n J_n(k rho) exp(i n phi)
    F_sca = sum_n i^n a_n H_n^(1)(k rho) exp(i n phi)
    F_int = sum_n i^n c_n J_n(m k rho) exp(i n phi)      (rho < R)

Matching tangential fields at ``rho = R`` (``F`` continuous; ``(1/eps)
dF/drho`` continuous for TE and ``dF/drho`` continuous for TM) gives::

    a_n = [ g J_n'(m x) J_n(x) - J_n(m x) J_n'(x) ]
          / [ J_n(m x) H_n'(x) - g J_n'(m x) H_n(x) ],   x = k R

with ``g = m`` (TM) or ``g = 1/m`` (TE), ``m = sqrt(eps_r)`` (complex allowed).
For a perfect electric conductor: ``a_n = -J_n(x)/H_n(x)`` (TM) and
``a_n = -J_n'(x)/H_n'(x)`` (TE).

Far-field observable
--------------------
The two-dimensional *echo width* (a.k.a. scattering width, the 2-D analogue
of radar cross section; units of length)::

    W(phi) = lim_{rho->inf} 2 pi rho |F_sca|^2 / |F_inc|^2
           = (4/k) | sum_n a_n exp(i n phi) |^2

The total scattering and extinction widths are::

    sigma_s = (1/2pi) \\oint W dphi = (4/k) sum_n |a_n|^2
    sigma_e = -(4/k) Re sum_n a_n        (optical theorem)

``phi`` is measured from the +x axis (the propagation direction), so
``phi = 0`` is forward scattering and ``phi = pi`` is backscatter.

Note: this is *not* the three-dimensional radar cross section (units of area).
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from scipy import special

Polarization = Literal["TM", "TE"]

__all__ = [
    "mie_coefficients",
    "echo_width",
    "total_scattering_width",
    "total_extinction_width",
    "scattering_efficiency",
    "backscatter_width",
    "near_field",
    "default_n_max",
]


def default_n_max(size_parameter: float, m: complex = 1.0) -> int:
    """Series truncation order (Wiscombe-type rule, padded for safety).

    Uses the larger of the exterior and interior size parameters.
    """
    x = abs(size_parameter) * max(1.0, abs(m))
    return int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 10.0))


def _validate(size_parameter: float, pol: str, material: str) -> None:
    if not np.isfinite(size_parameter) or size_parameter <= 0:
        raise ValueError(f"size_parameter must be positive and finite, got {size_parameter!r}")
    if pol not in ("TM", "TE"):
        raise ValueError(f"polarization must be 'TM' or 'TE', got {pol!r}")
    if material not in ("dielectric", "pec"):
        raise ValueError(f"material must be 'dielectric' or 'pec', got {material!r}")


def mie_coefficients(
    size_parameter: float,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
    n_max: int | None = None,
) -> np.ndarray:
    """Scattering coefficients ``a_n`` for ``n = 0 .. n_max``.

    Because ``a_{-n} = a_n`` only non-negative orders are returned; use
    :func:`echo_width` for sums over the full two-sided series.

    Parameters
    ----------
    size_parameter : ``x = k R`` with ``k`` the vacuum wavenumber.
    eps_r : relative permittivity (may be complex for lossy media); ignored
        for ``material="pec"``.
    polarization : ``"TM"`` (E_z) or ``"TE"`` (H_z).
    material : ``"dielectric"`` or ``"pec"``.
    n_max : highest order retained (default :func:`default_n_max`).
    """
    _validate(size_parameter, polarization, material)
    x = float(size_parameter)
    m = np.sqrt(complex(eps_r)) if material == "dielectric" else 1.0
    if material == "dielectric" and eps_r == 0:
        raise ValueError("eps_r must be non-zero")
    if n_max is None:
        n_max = default_n_max(x, m)
    n = np.arange(0, n_max + 1)

    jx, jpx = special.jv(n, x), special.jvp(n, x)
    hx, hpx = special.hankel1(n, x), special.h1vp(n, x)

    if material == "pec":
        if polarization == "TM":
            return -jx / hx
        return -jpx / hpx

    mx = m * x
    jm, jpm = special.jv(n, mx), special.jvp(n, mx)
    g = m if polarization == "TM" else 1.0 / m
    num = g * jpm * jx - jm * jpx
    den = jm * hpx - g * jpm * hx
    return num / den


def _full_series_sum(a: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """``sum_{n=-N}^{N} a_|n| exp(i n phi)`` using the symmetry ``a_-n = a_n``."""
    phi = np.asarray(phi, dtype=float)
    n = np.arange(1, a.size)
    out = a[0] + 2.0 * np.tensordot(np.cos(np.multiply.outer(phi, n)), a[1:], axes=([-1], [0]))
    return out


def echo_width(
    phi: np.ndarray | float,
    size_parameter: float,
    radius: float = 1.0,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
    n_max: int | None = None,
) -> np.ndarray:
    """2-D echo width ``W(phi) = (4/k)|sum a_n e^{i n phi}|^2`` (units of ``radius``).

    ``k = size_parameter / radius``.  Returns an array shaped like ``phi``.
    """
    a = mie_coefficients(size_parameter, eps_r, polarization, material, n_max)
    k = size_parameter / radius
    s = _full_series_sum(a, np.atleast_1d(phi))
    w = (4.0 / k) * np.abs(s) ** 2
    return w if np.ndim(phi) else float(w[0])


def total_scattering_width(
    size_parameter: float,
    radius: float = 1.0,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
    n_max: int | None = None,
) -> float:
    """``sigma_s = (4/k) sum_n |a_n|^2`` (length)."""
    a = mie_coefficients(size_parameter, eps_r, polarization, material, n_max)
    k = size_parameter / radius
    return float((4.0 / k) * (abs(a[0]) ** 2 + 2.0 * np.sum(np.abs(a[1:]) ** 2)))


def total_extinction_width(
    size_parameter: float,
    radius: float = 1.0,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
    n_max: int | None = None,
) -> float:
    """``sigma_e = -(4/k) Re sum_n a_n`` (length); equals ``sigma_s`` if lossless."""
    a = mie_coefficients(size_parameter, eps_r, polarization, material, n_max)
    k = size_parameter / radius
    return float(-(4.0 / k) * np.real(a[0] + 2.0 * np.sum(a[1:])))


def scattering_efficiency(
    size_parameter: float,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
) -> float:
    """``Q_sca = sigma_s / (2 R)`` (dimensionless geometric-width efficiency)."""
    return total_scattering_width(size_parameter, 1.0, eps_r, polarization, material) / 2.0


def backscatter_width(
    size_parameter: float,
    radius: float = 1.0,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
) -> float:
    """Monostatic (``phi = pi``) echo width."""
    return float(echo_width(np.pi, size_parameter, radius, eps_r, polarization, material))


def near_field(
    x: np.ndarray,
    y: np.ndarray,
    size_parameter: float,
    radius: float = 1.0,
    eps_r: complex = 4.0,
    polarization: Polarization = "TM",
    material: str = "dielectric",
    n_max: int | None = None,
) -> np.ndarray:
    """Complex total axial field ``F(x, y)`` (incident amplitude 1 at the origin).

    Inside the cylinder the interior expansion is used; outside, incident plus
    scattered.  For ``material="pec"`` the interior field is returned as 0.
    ``x`` and ``y`` are in the same units as ``radius``.
    """
    _validate(size_parameter, polarization, material)
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    rho = np.hypot(xs, ys)
    phi = np.arctan2(ys, xs)
    k = size_parameter / radius
    m = np.sqrt(complex(eps_r)) if material == "dielectric" else 1.0
    if n_max is None:
        n_max = default_n_max(size_parameter, m)
    n = np.arange(0, n_max + 1)
    a = mie_coefficients(size_parameter, eps_r, polarization, material, n_max)
    x0 = size_parameter

    field = np.zeros(rho.shape, dtype=complex)
    outside = rho >= radius
    # orders n>=1 appear twice (+n and -n): i^n e^{i n phi} + i^n e^{-i n phi} * (-1)^n/(-1)^n
    # J_{-n} = (-1)^n J_n,  H_{-n} = (-1)^n H_n  =>  the -n term equals the +n term with
    # e^{-i n phi}; summed they give 2 cos(n phi).
    fac = np.where(n == 0, 1.0, 2.0) * (1j) ** n

    if np.any(outside):
        r = rho[outside]
        ph = phi[outside]
        jn = special.jv(n[:, None], k * r[None, :])
        hn = special.hankel1(n[:, None], k * r[None, :])
        cosn = np.cos(n[:, None] * ph[None, :])
        field[outside] = np.sum(fac[:, None] * (jn + a[:, None] * hn) * cosn, axis=0)

    inside = ~outside
    if np.any(inside) and material == "dielectric":
        r = rho[inside]
        ph = phi[inside]
        jx, jpx = special.jv(n, x0), special.jvp(n, x0)
        hx, hpx = special.hankel1(n, x0), special.h1vp(n, x0)
        mx = m * x0
        jm, jpm = special.jv(n, mx), special.jvp(n, mx)
        c = (jx + a * hx) / jm
        jn_in = special.jv(n[:, None], m * k * r[None, :])
        cosn = np.cos(n[:, None] * ph[None, :])
        field[inside] = np.sum(fac[:, None] * c[:, None] * jn_in * cosn, axis=0)
        del jpx, hpx, jpm

    return field
