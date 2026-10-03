# Future work

The present repository covers verified two-dimensional electromagnetic scattering from a dielectric cylinder and resonance calculations for a one-dimensional periodic cylinder array. The items below are natural extensions of that work; they are not implemented unless noted otherwise in the current source tree.

## 1. Electromagnetic extensions

- Resolve the high-$Q$ TM feature near $\Lambda\approx2.4$–2.6 with a dedicated convergence study, finer parameter sampling, and field-symmetry analysis.
- Test the expected near-BIC scaling of $Q$ with an appropriate tuning parameter, including a dedicated near-$\Gamma$ Bloch-wavenumber study where applicable.
- Compare the periodic-array spectra with an independent method such as rigorous coupled-wave analysis (RCWA) or lattice-sum / multiple-scattering theory.
- Complete the large-pitch limit and compare the lowest array resonance with the isolated-cylinder Mie resonance.
- Extend the geometry to two-dimensional crossed lattices and, later, full three-dimensional unit cells.
- Add material dispersion and loss.
- Study finite arrays to quantify edge and finite-aperture effects.
- Explore inverse design of the unit cell using adjoint optimization, with finite-difference gradient checks.

## 2. RF sensing and signal-processing extensions

The electromagnetic solver can also be used as the front end of a synthetic RF sensing chain. Possible additions include:

- broadband scattering sweeps beyond the present benchmark band;
- synthetic complex I/Q returns generated from scattered fields;
- transmitted waveforms and matched filtering;
- range estimation and range profiles;
- Doppler and range-Doppler processing;
- noise and clutter models;
- CFAR detection;
- antenna-array models;
- conventional beamforming;
- MVDR / adaptive nulling;
- an integrated electromagnetic-to-detection example.

## 3. Current boundary of the repository

The code currently stops at electromagnetic scattering, periodic-array resonance analysis, and the associated verification studies. The RF sensing chain, inverse design, independent RCWA cross-check, material loss/dispersion, and higher-dimensional lattices remain future work.


## References

Method references for the present implementation are collected in [`07_references.md`](07_references.md).
