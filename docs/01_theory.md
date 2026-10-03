# Theory

This document summarizes the analytical model used for the single-cylinder scattering calculations. Theory for the periodic cylinder array is described in [`04_periodic_structures.md`](04_periodic_structures.md).

## 1. Maxwell's equations and constitutive relations

In a linear, isotropic, non-magnetic, lossless dielectric ($\mu=\mu_0$),

$$
\nabla\times\mathbf E=-\mu_0\,\partial_t\mathbf H,\qquad
\nabla\times\mathbf H=\varepsilon_0\varepsilon_r(\mathbf r)\,\partial_t\mathbf E,\qquad
\mathbf D=\varepsilon_0\varepsilon_r\mathbf E,\ \ \mathbf B=\mu_0\mathbf H .
$$

MEEP uses scale-free units with $c=\varepsilon_0=\mu_0=1$. A length $L$ is arbitrary, frequency is $f=\omega/2\pi$ in units of $c/L$, and the vacuum wavenumber is $k=2\pi f$. Results are therefore expressed in terms of the size parameter $ka=2\pi f a/c$, where $a$ is the cylinder radius. The same dimensionless calculation can be mapped to any physical scale; for example, $a=1$ mm gives $f=ka\times47.7$ GHz.

## 2. Two-dimensional reduction

For a structure invariant along $z$ and illuminated at normal incidence, with the wave vector in the $xy$ plane, Maxwell's equations separate into two scalar problems. Let $F$ denote the axial field component:

| Name in this repo | Axial field $F$ | Also called | Interface conditions at $\rho=a$ |
|---|---|---|---|
| `"TM"` | $E_z$ | $E$-polarization, *s* | $F$ and $\partial_\rho F$ continuous |
| `"TE"` | $H_z$ | $H$-polarization, *p* | $F$ and $\varepsilon_r^{-1}\partial_\rho F$ continuous |

Away from sources, both satisfy

$$
(\nabla_\perp^2+k^2\varepsilon_r)F=0.
$$

> **Naming note.** In this repository, TM/TE are defined relative to the cylinder axis ($z$), following the usual convention for scattering by cylinders: TM$_z$ corresponds to $E_z$ and TE$_z$ to $H_z$.

## 3. Exact solution for a circular cylinder

The cylindrical-wave treatment follows the classical dielectric-cylinder scattering formulation of [Wait (1955), *Scattering of a Plane Wave from a Circular Dielectric Cylinder at Oblique Incidence*](https://doi.org/10.1139/p55-024), specialized here to normal incidence. Take a unit-amplitude plane wave travelling along $+x$ with time dependence $e^{-i\omega t}$. The incident, scattered, and interior fields are expanded as

$$
F_\text{inc}=e^{ikx}=\sum_{n=-\infty}^{\infty} i^n J_n(k\rho)e^{in\phi},\quad
F_\text{sca}=\sum_n i^n a_n H^{(1)}_n(k\rho)e^{in\phi},\quad
F_\text{int}=\sum_n i^n c_n J_n(mk\rho)e^{in\phi},
$$

with $m=\sqrt{\varepsilon_r}$. Matching the interface conditions at $\rho=a$ gives, with $x=ka$,

$$
a_n=\frac{g\,J_n'(mx)J_n(x)-J_n(mx)J_n'(x)}{J_n(mx)H_n^{(1)\prime}(x)-g\,J_n'(mx)H^{(1)}_n(x)},\qquad
g=\begin{cases} m & \text{TM}\\ 1/m&\text{TE}.\end{cases}
$$

For a perfect conductor,

$$
a_n=-\frac{J_n(x)}{H_n^{(1)}(x)}\qquad\text{(TM)},
$$

or

$$
a_n=-\frac{J_n'(x)}{H_n^{(1)\prime}(x)}\qquad\text{(TE)}.
$$

Since $J_{-n}=(-1)^nJ_n$ and the same relation holds for $H_n^{(1)}$, the coefficients satisfy $a_{-n}=a_n$. The implementation is in [`em_rf_lab.analytical`](../src/em_rf_lab/analytical.py).

## 4. Two-dimensional scattering width

In two dimensions the far field decays as $\rho^{-1/2}$, so the natural scattering observable has units of length. The echo width, also called the scattering width, is

$$
W(\phi)=\lim_{\rho\to\infty}2\pi\rho\,\frac{|F_\text{sca}|^2}{|F_\text{inc}|^2}
=\frac{4}{k}\Big|\sum_{n=-\infty}^{\infty}a_n e^{in\phi}\Big|^2 .
$$

The total scattering and extinction widths are

$$
\sigma_s=\frac1{2\pi}\oint W\,d\phi=\frac4k\sum_n|a_n|^2,\qquad
\sigma_e=-\frac4k\,\mathrm{Re}\sum_n a_n .
$$

These are two-dimensional quantities. A three-dimensional radar cross section has units of area and is not the same observable. Here $\phi=0$ is forward scattering and $\phi=\pi$ is backscatter.

Checks used in the analytical tests include:

- **Optical theorem / energy conservation:** $\sigma_e=\sigma_s$ for lossless scatterers, verified to $10^{-12}$ in [`tests/unit/test_analytical.py`](../tests/unit/test_analytical.py); $\sigma_e>\sigma_s$ when $\varepsilon_r$ is complex.
- **Rayleigh limit** for a TM dielectric cylinder: $\sigma_s\to\frac{\pi^2}{4}(\varepsilon_r-1)^2k^3a^4$.
- **Geometric-optics limit** for a PEC cylinder: $\sigma_s/(2a)\to2$.
- Direct numerical solution of the $2\times2$ boundary-matching system reproduces the closed-form coefficients.
- The interface conditions are checked by finite differences applied to the evaluated near field.

## 5. Dimensionless scaling

For the circular-cylinder problem, the solution depends only on $ka$, $\varepsilon_r$, and polarization. The normalized quantities $W/a$ and $\sigma_s/a$ are therefore universal functions of those parameters.

The calculation used here covers $ka\in[2,6]$ for $\varepsilon_r=4$ ($m=2$), corresponding to interior size parameters $mka\in[4,12]$. Several narrow resonances occur toward the upper end of this interval, including whispering-gallery-like modes that make the numerical comparison more demanding.
