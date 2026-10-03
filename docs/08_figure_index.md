# Figure index

This index collects every figure used in the documentation in numerical order. Figure numbering is continuous across documents rather than restarting in each file. The images below are the same PNG files embedded in the numbered documentation files.

## Fig. 1 — Computational domain and monitor geometry

![Fig. 1 — Computational domain and monitor geometry](../figures/tier1_geometry.png)

**Appears in:** `02_numerical_method.md`, under **Setup of the cylinder-scattering calculation**, immediately after the configuration table.

Computational domain, PML, source line, dielectric cylinder, and closed near-to-far/flux box.

## Fig. 2 — FDTD and exact Mie angular scattering patterns

![Fig. 2 — FDTD and exact Mie angular scattering patterns](../figures/tier1_patterns.png)

**Appears in:** `03_verification_validation.md`, under **1. Cylinder scattering versus exact Mie theory**, before §1.1.

FDTD and exact Mie angular echo-width patterns at four size parameters.

## Fig. 3 — Grid convergence, numerical dispersion, and self-convergence

![Fig. 3 — Grid convergence, numerical dispersion, and self-convergence](../figures/tier1_convergence.png)

**Appears in:** `03_verification_validation.md`, under **§1.1 Grid convergence**, after the convergence table and discussion.

Grid convergence, fitted spectral shift due to numerical dispersion, and self-convergence between successive grids.

## Fig. 4 — Numerical-setup sensitivity

![Fig. 4 — Numerical-setup sensitivity](../figures/tier1_domain_studies.png)

**Appears in:** `03_verification_validation.md`, under **§1.2 Domain, source, and duration sensitivity**, after the sensitivity table.

Sensitivity to PML thickness, Gaussian beam waist, near-to-far box size, run time, and source bandwidth.

## Fig. 5 — Accuracy-cost tradeoff

![Fig. 5 — Accuracy-cost tradeoff](../figures/tier1_error_vs_runtime.png)

**Appears in:** `03_verification_validation.md`, under **§1.3 Run-duration artifact**.

Accuracy versus wall-clock cost for the resolution ladder.

## Fig. 6 — Total scattering-width spectrum and signed error

![Fig. 6 — Total scattering-width spectrum and signed error](../figures/tier1_spectrum.png)

**Appears in:** `03_verification_validation.md`, under **§1.4 Additional checks**.

Total scattering-width spectrum from two independent FDTD estimators compared with exact Mie theory, together with signed relative error.

## Fig. 7 — Near-field comparison with the exact Mie solution

![Fig. 7 — Near-field comparison with the exact Mie solution](../figures/tier1_field_map.png)

**Appears in:** `03_verification_validation.md`, under **§1.4 Additional checks**, after Fig. 6.

FDTD field, exact Mie field, and pointwise absolute difference near the dielectric cylinder.

## Fig. 8 — Periodic-cell verification against the exact slab solution

![Fig. 8 — Periodic-cell verification against the exact slab solution](../figures/tier2_slab_validation.png)

**Appears in:** `03_verification_validation.md`, under **2. Periodic-cell verification with an exact slab**.

Exact-versus-FDTD slab spectra, grid convergence, and oblique-incidence residuals.

## Fig. 9 — Sensitivity of periodic-array resonance metrics

![Fig. 9 — Sensitivity of periodic-array resonance metrics](../figures/tier2_convergence.png)

**Appears in:** `03_verification_validation.md`, under **3. Convergence of periodic-array resonance metrics**.

Sensitivity of resonance frequency and quality factor to resolution and other numerical setup choices.

## Fig. 10 — Normal-incidence reflectance versus pitch

![Fig. 10 — Normal-incidence reflectance versus pitch](../figures/tier2_pitch_maps.png)

**Appears in:** `05_results.md`, under **§1.1 Coarse reflectance survey**.

TM and TE normal-incidence reflectance maps versus lattice pitch and frequency.

## Fig. 11 — Lowest strong reflection resonance versus lattice pitch

![Fig. 11 — Lowest strong reflection resonance versus lattice pitch](../figures/tier2_lowest_resonance.png)

**Appears in:** `05_results.md`, under **§1.2 Lowest strong reflection resonance**.

Resonance frequency, quality factor, and usable reflection bandwidth versus pitch.

## Fig. 12 — TM reflectance and Harminv ringdown modes versus Bloch wavenumber

![Fig. 12 — TM reflectance and Harminv ringdown modes versus Bloch wavenumber](../figures/tier2_bands_TM_L3.png)

**Appears in:** `05_results.md`, under **3. Bloch-wavenumber dependence**.

TM reflectance map and Harminv ringdown modes versus normalized Bloch wavenumber at pitch 3.

## Numbering convention

The figure numbers above are the canonical documentation numbers. They are assigned in first-appearance order across the numbered Markdown documents:

1. `02_numerical_method.md`
2. `03_verification_validation.md`
3. `05_results.md`

Files without figures do not consume figure numbers. If a new figure is inserted earlier in the reading order, renumber subsequent figures and update this index at the same time.
