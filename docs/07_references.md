# References and sources

This bibliography gives the full publication details for sources cited inline throughout the documentation. The original build transcript did **not** maintain a contemporaneous bibliography, so these references should be read as the canonical literature supporting the methods used in the project, not as a claim that every paper was consulted line-by-line during implementation.

## Core references

1. **J. R. Wait**, “Scattering of a Plane Wave from a Circular Dielectric Cylinder at Oblique Incidence,” *Canadian Journal of Physics*, **33**(5), 189–195 (1955). DOI: 10.1139/p55-024.
   Supports the exact cylindrical-wave treatment of electromagnetic scattering by an infinite dielectric circular cylinder. The present benchmark uses the normal-incidence specialization.

2. **K. S. Yee**, “Numerical Solution of Initial Boundary Value Problems Involving Maxwell's Equations in Isotropic Media,” *IEEE Transactions on Antennas and Propagation*, **14**(3), 302–307 (1966). DOI: 10.1109/TAP.1966.1138693.
   Original staggered-grid finite-difference time-domain formulation used as the basis of the Yee scheme.

3. **A. Farjadpour, D. Roundy, A. Rodriguez, M. Ibanescu, P. Bermel, J. D. Joannopoulos, S. G. Johnson, and G. W. Burr**, “Improving Accuracy by Subpixel Smoothing in the Finite-Difference Time Domain,” *Optics Letters*, **31**(20), 2972–2974 (2006). DOI: 10.1364/OL.31.002972.
   Describes subpixel smoothing for discontinuous dielectric interfaces and its effect on FDTD convergence.

4. **A. F. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J. D. Joannopoulos, and S. G. Johnson**, “Meep: A Flexible Free-Software Package for Electromagnetic Simulations by the FDTD Method,” *Computer Physics Communications*, **181**(3), 687–702 (2010). DOI: 10.1016/j.cpc.2009.11.008.
   Primary software-method reference for MEEP, including FDTD implementation details, subpixel treatment, scripting, and signal processing.

5. **V. A. Mandelshtam and H. S. Taylor**, “Spectral Analysis of Time Correlation Function for a Dissipative Dynamical System Using Filter Diagonalization: Application to Calculation of Unimolecular Decay Rates,” *Physical Review Letters*, **78**, 3274–3277 (1997). DOI: 10.1103/PhysRevLett.78.3274.
   Methodological basis for harmonic inversion / filter diagonalization used by Harminv to extract resonant frequencies and decay rates from time-domain signals.

6. **S. Fan and J. D. Joannopoulos**, “Analysis of Guided Resonances in Photonic Crystal Slabs,” *Physical Review B*, **65**, 235112 (2002). DOI: 10.1103/PhysRevB.65.235112.
   Background for leaky guided resonances, dispersion, symmetry effects, and resonance quality factors in periodic dielectric structures.

7. **S. Fan, W. Suh, and J. D. Joannopoulos**, “Temporal Coupled-Mode Theory for the Fano Resonance in Optical Resonators,” *Journal of the Optical Society of America A*, **20**(3), 569–572 (2003). DOI: 10.1364/JOSAA.20.000569.
   Basis for the single-resonance coupled-mode/Fano line-shape model used to fit periodic-array reflection features.

8. **C. W. Hsu, B. Zhen, A. D. Stone, J. D. Joannopoulos, and M. Soljačić**, “Bound States in the Continuum,” *Nature Reviews Materials*, **1**, 16048 (2016). DOI: 10.1038/natrevmats.2016.48.
   Review reference for symmetry-protected and accidental bound states in the continuum and their radiative-decoupling mechanisms.

9. **A. Hessel and A. A. Oliner**, “A New Theory of Wood's Anomalies on Optical Gratings,” *Applied Optics*, **4**(10), 1275–1297 (1965). DOI: 10.1364/AO.4.001275.
   Classical reference distinguishing Rayleigh-order anomalies from resonance anomalies associated with guided or surface-wave coupling in periodic gratings.

## Software and documentation

- **MEEP documentation:** https://meep.readthedocs.io/
  Used for solver-specific APIs such as near-to-far transforms, DFT monitors, Bloch-periodic boundaries, flux normalization, symmetries, and Harminv integration.

## How these references map to the repository

| Repository topic | Main references |
|---|---|
| Exact infinite-cylinder scattering | [1] |
| Yee-grid FDTD and grid convergence | [2] |
| Dielectric-interface subpixel smoothing | [3], [4] |
| MEEP implementation and numerical workflow | [4] |
| Harminv / harmonic inversion | [5] |
| Periodic guided/lattice resonances | [6] |
| Coupled-mode and Fano line-shape fitting | [7] |
| Bound states in the continuum | [8] |
| Rayleigh-Wood anomalies | [9] |
