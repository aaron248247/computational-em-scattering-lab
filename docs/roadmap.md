# Roadmap

Status legend: **COMPLETE** = implemented, tested and with results in this repository; **IN PROGRESS** = partly
implemented; **PLANNED** = not started (design intent only). Nothing marked PLANNED exists in the code.

## Delivered in this release

| ID | Item | Status |
|---|---|---|
| V1 | Exact cylinder (Mie) reference solution, TM/TE, dielectric and PEC | **COMPLETE** |
| V2 | FDTD scattering workflow: incident normalisation, scattered-field subtraction, near-to-far transform | **COMPLETE** |
| V3 | Verification against exact theory: grid convergence, numerical-dispersion shift, PML / domain / finite-beam / run-time sensitivity | **COMPLETE** |
| V4 | Mirror-symmetry reduction verified against the full domain (TM and TE) | **COMPLETE** |
| P1 | Periodic-layer reflectance/transmittance workflow (Bloch boundaries, oblique incidence), verified against the exact slab solution | **COMPLETE** |
| P2 | Pitch and polarisation study of a dielectric-cylinder array (resonance frequency, Q, bandwidth) | **COMPLETE** (see `docs/results.md`) |
| P3 | Angle / Bloch-wavenumber dependence (R/T maps + ringdown modes) | see `docs/results.md` |
| P4 | Resolution and domain convergence of the resonance metrics | see `docs/results.md` |

## Next: toward RF sensing and signal processing (all PLANNED)

The original sensing roadmap is retained. The two tiers above are the electromagnetic foundation
(validated scattering and resonant-layer response) on which these build.

| ID | Milestone | Status |
|---|---|---|
| M2 | Broadband scattering and frequency sweeps (partly covered: the benchmark already yields 41-point spectra) | PLANNED (extension) |
| M3 | Synthetic complex I/Q returns from scattered fields | PLANNED |
| M4 | Transmitted waveforms and matched filtering | PLANNED |
| M5 | Range estimation / range profiles | PLANNED |
| M6 | Doppler and range–Doppler processing | PLANNED |
| M7 | Noise and clutter models | PLANNED |
| M8 | CFAR detection | PLANNED |
| M9 | Antenna-array model | PLANNED |
| M10 | Conventional beamforming | PLANNED |
| M11 | MVDR / adaptive nulling | PLANNED |
| M12 | Integrated EM-to-detection experiment | PLANNED |

## Electromagnetic extensions of the primary investigation (all PLANNED)

* Two-dimensional (crossed) lattices and full 3-D unit cells; material dispersion and loss.
* Inverse design of the unit cell (adjoint optimisation, finite-difference gradient checks).
* Comparison with rigorous coupled-wave analysis (RCWA) or lattice-sum (multiple-scattering) theory as a
  second independent reference for the array spectra.
* Finite-aperture arrays and edge effects.
