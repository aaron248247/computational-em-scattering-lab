# Numerical method

All simulations use [MEEP](https://meep.readthedocs.io) (finite-difference time-domain, FDTD), driven by the
code in `src/em_rf_lab`. This document describes the method; the evidence that it is accurate is in
[`verification_validation.md`](verification_validation.md).

## FDTD in one paragraph

Maxwell's curl equations are discretised on a Yee grid: $\mathbf E$ and $\mathbf H$ components live on
staggered positions and are advanced leap-frog in time with step $\Delta t=S\,\Delta x$ ($c=1$), where MEEP's
default Courant number is $S=0.5$ (the 2-D stability limit is $1/\sqrt2$). The scheme is second-order accurate
in space and time for smooth fields, so the phase-velocity error (numerical dispersion) of a wave with $N$
points per wavelength scales as $(\pi/N)^2/6$ or better, which is why resolution must be set by the shortest
*in-medium* wavelength: here $\lambda_\min/n=2\pi/(n\,k_{\max})=1/(2f_\max)\approx0.52$ length units, i.e.
$\approx10$ points per wavelength at resolution 20.

Material interfaces are treated with MEEP's **subpixel smoothing** (anisotropic effective permittivity at
boundary pixels), which restores second-order convergence for a curved dielectric boundary. This is
deliberately why a *dielectric* cylinder, not a PEC (staircased, first-order), is used as the benchmark.

## Setup of the benchmark (Tier 1)

| Item | Choice |
|---|---|
| Dimensionality | 2-D ($xy$), structure invariant in $z$ |
| Scatterer | circular cylinder, $a=1$, $\varepsilon_r=4$ (lossless) |
| Polarisation | TM ($E_z$, primary) and TE ($H_z$) |
| Source | broadband Gaussian-pulse **line source** at $x=-4.5$, spanning the full cell height; amplitude profile $\exp(-(y/w)^2)$, $w=20$ (finite-beam approximation of a plane wave) |
| Frequencies | $ka\in[2,6]$, 41 uniformly spaced DFT frequencies |
| Boundaries | PML, thickness 5, all four sides |
| Cell | $20\times56$ (including PML) |
| Near-to-far surface | closed square box of half-width 2.5 around the cylinder (also the scattered-flux box) |
| Observation | far field at $\rho=10^5$ (evaluated with MEEP's exact 2-D Green's function), 91 angles on $[0,\pi]$ |
| Termination | stop when the axial field at an interior probe has decayed by `decay_tol` ($10^{-5}$) relative to its peak |
| Symmetry | mirror $y\to-y$ (phase $+1$ for $E_z$, $-1$ for $H_z$) — verified to reproduce the full-domain run to $10^{-13}$ |

## Workflow and normalisation

1. **Reference run** (vacuum). Records (i) the complex incident amplitude $F_\text{inc}(\omega)$ at the origin
   (DFT), (ii) the incident near-field and flux data on the closed monitor box, (iii) the incident profile
   along $y$ at $x=0$ (plane-wave quality diagnostic).
2. **Scattering run** (cylinder). The stored incident data are **subtracted** (`load_minus_near2far_data`,
   `load_minus_flux_data`), so the monitors contain only the scattered field.
3. **Echo width** from the near-to-far transform of the scattered near field:
   $W=2\pi\rho\,|F_s(\rho,\phi)|^2/|F_\text{inc}(0)|^2$.
4. **Total scattering width** by an independent route: scattered Poynting flux through the box divided by the
   incident intensity $I_0=\tfrac12|F_\text{inc}|^2$. MEEP's DFT flux is $\mathrm{Re}\,(\mathbf E\times\mathbf H^*)$
   with no factor $1/2$, hence $\sigma_s^\text{flux}=P_s/(2I_0)$; this convention is verified in every run by
   comparing the flux through a transverse line to the integral of the measured incident intensity profile
   (`flux_norm_check`, which converges to 1.000).

The agreement of the two estimators ($\sigma_s$ from the far-field integral vs from the scattered flux) is a
built-in internal-consistency check of the near-to-far transform.

## Quasi-plane-wave incidence and its error

A finite line source cannot produce an infinite plane wave. A hard-edged source suffers aperture-edge
diffraction (Fresnel number $N_F=(D/2)^2/\lambda z$); a Gaussian-tapered source removes that but leaves an
$O((a/w)^2)$ amplitude variation across the scatterer and finite angular spread $O(1/kw)$. Both are
quantified in the verification studies (beam waist $w$ and cell size) rather than assumed negligible.

## Run-time / spectral resolution

A DFT accumulated over a finite time $T$ has frequency resolution $\sim1/T$. Isolated high-$Q$ resonances
(quality factor $Q$) ring for $\sim Q/\pi f$ time units, so narrow spectral features require longer runs than
broad ones, and the decay test must probe a location where the slowest mode has appreciable field (see the
run-time study in the verification document).

## Reproducibility

Every run is identified by a hash of its complete configuration
([`SimulationConfig`](../src/em_rf_lab/config.py)); the run directory stores `config.json`, `meta.json`
(timestamps, git commit, MEEP/NumPy/SciPy versions, wall time, termination status) and `arrays.npz` (all
figure data). Raw MEEP HDF5 output is never stored.
