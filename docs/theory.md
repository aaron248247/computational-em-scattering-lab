# Theory

This document states the physics and the exact reference solution used in the **verification benchmark**
(Tier 1). The periodic-structure theory used by the primary investigation (Tier 2) is in
[`periodic_structures.md`](periodic_structures.md).

## 1. Maxwell's equations and constitutive relations

In a linear, isotropic, non-magnetic, lossless dielectric ($\mu=\mu_0$):

$$
\nabla\times\mathbf E=-\mu_0\,\partial_t\mathbf H,\qquad
\nabla\times\mathbf H=\varepsilon_0\varepsilon_r(\mathbf r)\,\partial_t\mathbf E,\qquad
\mathbf D=\varepsilon_0\varepsilon_r\mathbf E,\ \ \mathbf B=\mu_0\mathbf H .
$$

MEEP works in scale-free units with $c=\varepsilon_0=\mu_0=1$: a length $L$ is arbitrary, frequency is
$f=\omega/2\pi$ in units of $c/L$, and the vacuum wavenumber is $k=2\pi f$. Results are therefore
stated in terms of the **size parameter** $ka=2\pi f a/c$ (with cylinder radius $a$) and are valid for any
physical scale (for example $a=1$ mm gives $f = ka \times 47.7$ GHz).

## 2. Two-dimensional reduction

For a structure that is invariant along $z$ and illuminated at normal incidence (wave vector in the
$xy$-plane) the fields split into two independent scalar problems. Let $F$ denote the axial component:

| Name in this repo | Axial field $F$ | Also called | Interface conditions at $\rho=a$ |
|---|---|---|---|
| `"TM"` | $E_z$ | $E$-polarisation, *s* | $F$ and $\partial_\rho F$ continuous |
| `"TE"` | $H_z$ | $H$-polarisation, *p* | $F$ and $\varepsilon_r^{-1}\partial_\rho F$ continuous |

Both satisfy the scalar Helmholtz equation $(\nabla_\perp^2+k^2\varepsilon_r)F=0$ away from sources.

> Naming caution: "TM/TE" here refers to fields *relative to the cylinder axis* ($z$), following the
> usual scattering-by-cylinders convention (TM$_z$ = $E_z$), not relative to a plane of incidence.

## 3. Exact solution: Mie series for a circular cylinder

Plane wave along $+x$, time dependence $e^{-i\omega t}$, unit amplitude:

$$
F_\text{inc}=e^{ikx}=\sum_{n=-\infty}^{\infty} i^n J_n(k\rho)e^{in\phi},\quad
F_\text{sca}=\sum_n i^n a_n H^{(1)}_n(k\rho)e^{in\phi},\quad
F_\text{int}=\sum_n i^n c_n J_n(mk\rho)e^{in\phi}
$$

with $m=\sqrt{\varepsilon_r}$. Matching the interface conditions gives, with $x=ka$,

$$
a_n=\frac{g\,J_n'(mx)J_n(x)-J_n(mx)J_n'(x)}{J_n(mx)H_n^{(1)\prime}(x)-g\,J_n'(mx)H^{(1)}_n(x)},\qquad
g=\begin{cases} m & \text{TM}\\ 1/m&\text{TE}\end{cases}
$$

and for a perfect conductor $a_n=-J_n(x)/H_n^{(1)}(x)$ (TM) or $-J_n'(x)/H_n^{(1)\prime}(x)$ (TE).
Since $J_{-n}=(-1)^nJ_n$ and likewise for $H^{(1)}$, $a_{-n}=a_n$.
Implemented in [`em_rf_lab.analytical`](../src/em_rf_lab/analytical.py).

## 4. Physical observable: the 2-D echo width (not "RCS")

In two dimensions the far field decays as $\rho^{-1/2}$, so the natural observable is a **length**, the
*echo width* (also "scattering width"):

$$
W(\phi)=\lim_{\rho\to\infty}2\pi\rho\,\frac{|F_\text{sca}|^2}{|F_\text{inc}|^2}
=\frac{4}{k}\Big|\sum_{n=-\infty}^{\infty}a_n e^{in\phi}\Big|^2 .
$$

The total scattering and extinction widths are

$$
\sigma_s=\frac1{2\pi}\oint W\,d\phi=\frac4k\sum_n|a_n|^2,\qquad
\sigma_e=-\frac4k\,\mathrm{Re}\sum_n a_n .
$$

The three-dimensional radar cross section has units of area and is a different quantity; this
repository never labels $W$ or $\sigma_s$ "RCS". $\phi=0$ is forward scattering and $\phi=\pi$ is
backscatter. Shape-independent checks used in the test suite:

* **Optical theorem / energy conservation:** $\sigma_e=\sigma_s$ for lossless scatterers (verified to $10^{-12}$ in
  [`tests/unit/test_analytical.py`](../tests/unit/test_analytical.py)); $\sigma_e>\sigma_s$ when $\varepsilon_r$ is complex.
* **Rayleigh limit** (TM dielectric): $\sigma_s\to\frac{\pi^2}{4}(\varepsilon_r-1)^2k^3a^4$.
* **Geometric-optics limit** (PEC): $\sigma_s/(2a)\to2$ (extinction paradox).
* Direct numerical solution of the $2\times2$ boundary-matching system reproduces the closed form.
* Interface conditions verified by finite differences on the evaluated near field.

## 5. Dimensionless scaling

Everything depends only on $ka$, $\varepsilon_r$ and polarisation; $W/a$ and $\sigma_s/a$ are universal
functions of those. The benchmark band $ka\in[2,6]$ for $\varepsilon_r=4$ ($m=2$) corresponds to interior
size parameters $mka\in[4,12]$, which contains several sharp resonances (whispering-gallery-like modes
near the upper end of the band) that make the test demanding.
