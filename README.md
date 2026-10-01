# computational-em-rf-lab

A small, verified computational-electromagnetics workflow built on the [MEEP](https://meep.readthedocs.io) FDTD solver:
a **graduate-level parametric study of the resonances of a periodic dielectric-cylinder array** (a one-dimensional
frequency-selective surface / metasurface), built on a **verification benchmark** in which the same FDTD workflow
reproduces the exact Mie solution for a single cylinder. Everything is generic, public and unclassified: canonical
dielectric cylinders in a dimensionless unit system (cylinder radius = 1, speed of light = 1).

| Part | Status |
|---|---|
| Primary investigation: resonances of a periodic cylinder array vs. pitch and polarisation | **COMPLETE** (core result); band-structure/ky scan **IN PROGRESS** (exploratory) |
| Verification benchmark: FDTD vs. exact Mie scattering of a dielectric cylinder | **COMPLETE** |
| Exact-slab verification of the periodic-cell workflow | **COMPLETE** |
| Radar signal-processing chain (I/Q, matched filtering, Doppler, CFAR, beamforming) | **PLANNED** — nothing of this exists in the code |
| Independent RCWA cross-check, 2-D lattices, inverse design | **PLANNED** |

## 1. Primary research investigation: collective resonances of a dielectric-cylinder array

**Question.** For an infinite row of lossless dielectric cylinders ($\varepsilon_r=4$, radius 1) with pitch $\Lambda$,
illuminated at normal incidence, how do the frequency $f_0$, quality factor $Q$ and usable (≥ 90 % reflection)
bandwidth of the lowest strong reflection resonance depend on $\Lambda$ and polarisation (TM $E_z$ / TE $H_z$)?
Is the response consistent with the physics (energy conservation, diffraction-order onset, coupled-mode
line shapes), and how well converged is it?

**Method.** One Bloch-periodic unit cell, plane-wave incidence, reflectance and transmittance from Poynting-flux
monitors normalised to an empty-cell run, line-shape fits with coupled-mode theory, and an independent
time-domain (Harminv) ringdown estimate of $f_0$ and $Q$. Theory: [`docs/periodic_structures.md`](docs/periodic_structures.md).

**Results** (details and exact numbers: [`docs/results.md`](docs/results.md)).

![lowest resonance vs pitch](figures/tier2_lowest_resonance.png)

* **TE** ($H_z$): a broad, **fully reflecting** band ($R\to1.000$, lossless) whose frequency falls from $f_0=0.285$ to
  $0.178$ as $\Lambda$ goes from 2.2 to 5. $Q$ is low ($\approx6$ – $28$) and the fractional ≥ 90 % reflection band has a
  maximum of ≈ 6 % near $\Lambda\approx2.8$.
* **TM** ($E_z$): $Q$ is strongly non-monotonic. The lowest TM resonance forms a smooth branch in $f_0$, but $Q$ rises
  from ≈ 4 at $\Lambda=5$ to ≈ 65 at $\Lambda=3$, then to **~10³** around $\Lambda=2.4$ – $2.6$ and falls again to ≈ 180 at
  $\Lambda=2.8$ — the signature of a mode that nearly decouples from the radiation continuum. We have *not* proven this is
  a bound state in the continuum; the $Q$ values at $\Lambda=2.4$/$2.6$ are approximate (the ringdown is only marginally
  converged).
* **Convergence.** $f_0$ is stable to $\le5\times10^{-4}$ relative and $Q$ to 0.2 % across grid resolution 20–60, PML
  thickness, structure-to-PML margin and monitor gap (two representative cases, `docs/verification_validation.md`).
  Line-shape and ringdown estimates of $Q$ agree to 7 % where both exist (TM, $\Lambda=3$: 65.0 vs 62.7).

![pitch maps](figures/tier2_pitch_maps.png)

**Limitations.** 2-D geometry, lossless non-dispersive dielectric, an infinite row (no edges); $Q\lesssim10$ values are
indicative because broad Fano-like lines are not described by a single Lorentzian; no second independent solver
(RCWA planned); the coarse survey in the pitch map under-samples the narrowest TM lines and some short-pitch TM runs
did not decay (flagged in `docs/results.md`).

## 2. Verification benchmark: FDTD scattering from a dielectric cylinder vs. exact Mie theory

Before using the workflow on a periodic structure, it was verified on the one geometry with an exact solution. A TM
plane wave (quasi-plane wave from a Gaussian-tapered line source) scatters from a cylinder ($ka\in[2,6]$); the
far-field **echo width** $W(\phi)$ (a 2-D quantity; deliberately *not* called an "RCS") is obtained by near-to-far
transformation of the scattered field and compared with the exact series in `src/em_rf_lab/analytical.py`.

![patterns](figures/tier1_patterns.png)

| Resolution (px/unit) | 10 | 15 | 20 | 30 | 40 | 60 |
|---|---|---|---|---|---|---|
| RMS relative $L^2$ error of $W(\phi)$ | 0.397 | 0.173 | 0.099 | 0.049 | 0.030 | 0.015 |

* Error decreases monotonically; fitted order 1.74. A rigid frequency shift from the Yee scheme's numerical
  dispersion explains most of it and converges with order 1.89.
* A first set of runs stalled at ≈ 1.4 % error; this was traced (run-time study) to **duration**, not discretisation —
  whispering-gallery modes near the rim ring long and are invisible to a centre probe. Both result sets are kept.
* PML thickness, beam waist, near-to-far box and source bandwidth are held to ≤ 0.5 % in the production setting.
* Half-domain mirror symmetry reproduces the full domain to $10^{-13}$.
* Worst-case error at resolution 60 is 6.1 % at the narrowest resonance ($ka\approx5.4$).

![convergence](figures/tier1_convergence.png)
![field map](figures/tier1_field_map.png)

## For technical reviewers

* The **question, physics, method, verification and limitations** of each part are in `docs/` (read
  `docs/verification_validation.md` first).
* Every figure is regenerated from stored results by `scripts/make_tier{1,2}_figures.py`; no figure is hand-made.
* Each run is stored under `results/` with `config.json`, `meta.json` (git commit, library versions, wall time,
  termination status) and `arrays.npz`; the directory name contains a hash of the full configuration.
* Disagreements and limitations are reported rather than tuned away: the duration artefact above, the
  under-resolved TM lines, near-grazing-incidence energy residuals in the slab check.
* The code is organised so the analytical reference (`analytical.py`), numerical workflow (`simulation.py`,
  `periodic.py`) and comparison (`validation.py`) are independent modules.

## Install and reproduce

```bash
# MEEP is not pip-installable; on Linux/macOS use conda-forge
conda create -n emrf -c conda-forge python=3.11 pymeep numpy scipy matplotlib pytest
conda activate emrf
pip install -e .

pytest -m "not meep"      # unit tests and stored-result checks, no MEEP needed (seconds)
pytest                    # everything incl. small MEEP integration tests (~4 min; MEEP tests auto-skip if MEEP is absent)
pytest -m results         # checks on the stored result files only

# Tier 1 (cylinder benchmark): ~30 min for resolutions 10-40 on 2 cores
python experiments/convergence/run_resolution_ladder.py --resolutions 10 15 20 30 40 60
python experiments/convergence/run_domain_studies.py
python experiments/baseline/run_baseline.py
python scripts/make_tier1_figures.py

# Tier 2 (periodic array)
python experiments/tier2/validate_slab.py
python experiments/tier2/pitch_study.py
python experiments/tier2/convergence_periodic.py
python experiments/tier2/resonance_vs_pitch.py --harminv
python experiments/tier2/highq_check.py
python experiments/tier2/band_structure.py --pitch 3.0       # exploratory
python scripts/make_tier2_figures.py
```

Tested with MEEP 1.34.0 (serial), Python 3.11, NumPy 2.4, SciPy 1.17, Matplotlib 3.10.
Wall times quoted in `results/*/meta.json`.

## Repository layout

```
src/em_rf_lab/   analytical, config, geometry, metrics, convergence, io, sources, monitors, farfield,
                 simulation, validation, plotting          (Tier 1)
                 periodic, periodic_sources, resonance, modes, tier2_analysis   (Tier 2)
experiments/     runnable studies (baseline, convergence, tier2)
scripts/         figure generation
tests/           unit (fast, no MEEP) and integration (MEEP; skipped if unavailable)
docs/            theory, numerical method, periodic structures, verification & validation, results, roadmap
results/         stored run data (small .npz/.json); figures/ - generated figures with captions
```

## Roadmap

Statuses in [`docs/roadmap.md`](docs/roadmap.md). Radar-style signal processing (synthetic I/Q, matched filtering,
range/Doppler, CFAR, array beamforming) is **PLANNED**; this repository currently contains none of it.

## Scope

Simple canonical geometries and dimensionless units only. No real platform, system, or sensitive parameter is
modelled.

## License and citation

MIT (see `LICENSE`). Citation metadata: `CITATION.cff`.
