# computational-em-scattering-lab

This repository contains a small computational-electromagnetics study built with the
[MEEP](https://meep.readthedocs.io) FDTD solver. The main calculation is a parameter
study of resonances in a one-dimensional periodic array of dielectric cylinders. The
same numerical workflow is first checked against the exact Mie-series solution for a
single dielectric cylinder. The geometry is deliberately simple: lossless dielectric cylinders with $\varepsilon_r=4$ in dimensionless units, with cylinder radius $R=1$ and $c=1$.

## Periodic dielectric-cylinder array

For an infinite row of cylinders with pitch $\Lambda$, the main quantities of interest
are the frequency $f_0$, quality factor $Q$, and the fractional bandwidth over which the
reflectance is at least 90%. Both TM ($E_z$) and TE ($H_z$) polarization are considered.

The periodic calculation uses a single Bloch-periodic unit cell with plane-wave
incidence. Reflectance and transmittance are obtained from Poynting-flux monitors and normalized to an empty-cell reference calculation. Resonances are characterized with line-shape fits and, where useful, with a separate time-domain Harminv ringdown estimate of $f_0$ and $Q$.

More background is in
[`docs/04_periodic_structures.md`](docs/04_periodic_structures.md), and the numerical results
are collected in [`docs/05_results.md`](docs/05_results.md).

![lowest resonance vs pitch](figures/tier2_lowest_resonance.png)

The TE branch is comparatively broad. Over the pitch range studied, the lowest strong
reflection resonance moves from about $f_0=0.285$ at $\Lambda=2.2$ to $f_0=0.178$ at
$\Lambda=5$. The corresponding $Q$ values are roughly 6--28, and the fractional
bandwidth with $R\ge0.9$ reaches about 6% near $\Lambda\approx2.8$.

The TM branch behaves differently. Its resonance frequency changes smoothly with pitch,
but its $Q$ is strongly non-monotonic: it is about 4 at $\Lambda=5$, about 65 at
$\Lambda=3$, and rises to order $10^3$ near $\Lambda=2.4$--$2.6$ before falling again.
This is consistent with a mode that becomes weakly coupled to the radiation continuum,
but the present calculations do not establish a bound state in the continuum. The
highest-$Q$ values at $\Lambda=2.4$ and $2.6$ come from long ringdown calculations and
are only marginally converged.

For two representative resonances, the fitted resonance frequency is stable to about
$5\times10^{-4}$ relative over the resolution and domain-parameter checks. The TM
$\Lambda=3$ case also provides a useful cross-check on $Q$: the fine line-shape fit gives
$Q\approx65.0$, while the Harminv ringdown estimate gives $Q\approx62.7$.

![pitch maps](figures/tier2_pitch_maps.png)

The periodic-array model is still limited to two dimensions, lossless nondispersive
materials, and an infinite row with no edge effects. Broad low-$Q$ Fano-like lines are
not always well represented by a single Lorentzian fit. The coarse pitch survey also
undersamples the narrowest TM resonances, which is why the short-pitch cases were
rechecked with finer frequency windows and longer ringdown runs. An independent RCWA
calculation has not yet been added.

## Single-cylinder verification

Before using the FDTD workflow for the periodic array, I checked it against the exact
solution for scattering from a single dielectric cylinder.

A TM quasi-plane wave from a Gaussian-tapered line source illuminates a cylinder over
$ka\in[2,6]$. The far-field echo width

$$
W(\phi)
$$

is obtained from MEEP's near-to-far transformation and compared with the exact Mie
series implemented in `src/em_rf_lab/analytical.py`. I use "echo width" rather than RCS because this is a two-dimensional scattering problem.

![patterns](figures/tier1_patterns.png)

| Resolution (px/unit)                  | 10    | 15    | 20    | 30    | 40    | 60    |
| ------------------------------------- | ----- | ----- | ----- | ----- | ----- | ----- |
| RMS relative $L^2$ error of $W(\phi)$ | 0.397 | 0.173 | 0.099 | 0.049 | 0.030 | 0.015 |

The error decreases monotonically with grid refinement. A fit over the finer grids gives
an observed convergence order of about 1.74. The fitted rigid spectral shift converges
with order about 1.89, close to the second-order behavior expected from numerical
dispersion in the Yee scheme.

The first version of the benchmark appeared to plateau near 1.4% error at high
resolution. That turned out to be a run-duration problem rather than a grid-resolution
problem: long-lived rim-localized modes were not being captured by the original decay
criterion. Moving the decay probe toward the cylinder rim and extending the minimum run time removed the plateau.

Separate checks were made for PML thickness, source waist, near-to-far-box placement, source bandwidth, and run duration. Half-domain mirror symmetry reproduces the full domain to numerical precision. At resolution 60, the largest single-frequency discrepancy is about 6.1% near the narrow resonance around $ka\approx5.4$.

![convergence](figures/tier1_convergence.png)

![field map](figures/tier1_field_map.png)

## Reproducibility

The plotting scripts read saved result files; they do not rerun the simulations. Runs
are stored under `results/` with the configuration, metadata, and numerical arrays needed to reproduce the comparisons. The analytical reference, FDTD workflow, and validation code are kept in separate modules.

Known discrepancies are kept in the results rather than tuned away. Examples include
the original short-run duration artifact, the under-resolved high-$Q$ TM lines in the
coarse pitch sweep, and the larger residuals near grazing incidence in the slab
validation.

## Install and run

```bash
# MEEP is not pip-installable; on Linux/macOS use conda-forge
conda create -n emrf -c conda-forge python=3.11 pymeep numpy scipy matplotlib pytest
conda activate emrf
pip install -e .

pytest -m "not meep"      # unit tests and stored-result checks, no MEEP needed
pytest                    # includes small MEEP integration tests
pytest -m results         # checks on stored result files only

# Single-cylinder benchmark
python experiments/convergence/run_resolution_ladder.py --resolutions 10 15 20 30 40 60
python experiments/convergence/run_domain_studies.py
python experiments/baseline/run_baseline.py
python scripts/make_tier1_figures.py

# Periodic array
python experiments/tier2/validate_slab.py
python experiments/tier2/pitch_study.py
python experiments/tier2/convergence_periodic.py
python experiments/tier2/resonance_vs_pitch.py --harminv
python experiments/tier2/highq_check.py
python experiments/tier2/band_structure.py --pitch 3.0
python scripts/make_tier2_figures.py
```

The saved results were generated with MEEP 1.34.0 (serial), Python 3.11,
NumPy 2.4, SciPy 1.17, and Matplotlib 3.10. Wall times are recorded in
`results/*/meta.json`.

## Repository layout

```text
src/em_rf_lab/   analytical, config, geometry, metrics, convergence, io, sources,
                 monitors, farfield, simulation, validation, plotting
                 periodic, periodic_sources, resonance, modes, tier2_analysis
experiments/     runnable studies for the cylinder benchmark and periodic array
scripts/         figure-generation scripts
tests/           unit and MEEP integration tests
docs/            theory, numerical method, periodic structures, verification,
                 results, and future work
results/         saved numerical runs
figures/         generated figures and figure notes
```

## Future work

Natural extensions are a more systematic study of the high-$Q$ TM branch, an independent
RCWA or lattice-sum cross-check, larger-pitch comparison with the isolated-cylinder
limit, two-dimensional or three-dimensional lattices, material loss and dispersion, and
inverse-design calculations.

A separate RF signal-processing extension—synthetic I/Q data, matched filtering,
range/Doppler processing, CFAR, and array processing—is not implemented in this
repository.

## Scope

The repository uses canonical geometries and dimensionless parameters only. It does not model a real platform, operational system, or sensitive parameter set.

## License and citation

MIT license; see `LICENSE`. Citation metadata are in `CITATION.cff`.
