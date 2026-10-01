# Verification and validation

*Verification* asks whether the equations are solved correctly; *validation* asks whether the right equations are
solved. Here both are tested against **exact analytical solutions** (Mie series for the cylinder, Airy/transfer
matrix for the slab), so the two overlap. All numbers below are produced by scripts in `experiments/` and stored in
`results/`; figures are produced by `scripts/make_tier*_figures.py`.

## Tier 1 — cylinder scattering vs. exact Mie theory (**COMPLETE**)

Configuration: TM, circular dielectric cylinder $\varepsilon_r=4$, $ka\in[2,6]$ (41 frequencies), 91 angles.
Metric: relative $L^2$ error of the echo width $W(\phi)$ per frequency, RMS over frequency.

### Grid convergence (time-converged runs)

| Resolution (px/unit) | RMS rel. L2 error of W | Error of $\sigma_s(ka)$ (L2) | Fitted dispersion shift $\delta$ |
|---|---|---|---|
| 10 | 0.397 | 0.201 | 2.25e-2 |
| 15 | 0.173 | 0.087 | 6.9e-3 |
| 20 | 0.099 | 0.054 | 4.2e-3 |
| 30 | 0.049 | 0.031 | 2.2e-3 |
| 40 | 0.030 | 0.021 | 1.4e-3 |
| 60 | 0.015 | 0.011 | 6.4e-4 |

* Error decreases monotonically; least-squares order for resolution ≥ 20 is **1.74** (local orders between successive resolutions: 2.05, 1.92, 1.75, 1.65, 1.72).
* The error is dominated by a *rigid fractional frequency shift* of the spectrum (numerical dispersion of the Yee
  scheme); the fitted shift converges with order **1.89** (≈ second order). Removing it lowers the RMS error by a
  factor of 1.5–2 depending on resolution — the remainder is the shape error.
* The maximum single-frequency error at resolution 60 is 6.1 %, at the narrowest resonance ($ka\approx5.4$), where
  a tiny frequency shift becomes a large amplitude error. Away from narrow resonances the error is far smaller
  (e.g. $ka=2$: 6e-4).
* Self-convergence (difference between successive grids) is reported in `figures/tier1_convergence.png`.

### Domain, source and duration sensitivity (resolution 30; differences relative to the production setting)

| Parameter | Values tested | Effect on echo-width error (RMS, vs production) |
|---|---|---|
| PML thickness | 2, 3, **5**, 8 | 1.3 %, 1.1 %, –, 0.3 % |
| Beam waist (finite-beam plane-wave error) | 5, 10, **20**, 30 | 4.4 %, 0.97 %, –, 0.34 % |
| Near-to-far box half-width | 2, **2.5**, 3.5, 5 | 0.4 %, –, 0.6 %, 0.8 % |
| Minimum run time | 50, 110, 250, **500**, 1000 | 5.3 %, 1.7 %, 0.39 %, 0.07 %, 0.01 % |
| Source bandwidth | 0.75, **1.5**, 3.0 × | 0.20 %, –, 0.01 % |

The production setting leaves these contributions ≤ 0.5 %, well below the 3–5 % discretisation error at
resolutions 30–40. (Restricted to $ka\lesssim3$ the sensitivity to a small beam waist or thin PML is larger — 4.9 % and 2.3 % RMS respectively — which is why the production
setting uses waist 20 and PML 5.)

### Duration artefact (an honest negative result kept in the repository)

A first ladder stopped when a field probe near the centre had decayed to $10^{-5}$ of its peak. Its error stopped
decreasing at ≈ 1.4 % for resolutions ≥ 60 (figure: grey squares). The cause is **not** discretisation: whispering-
gallery-type resonances near the rim ring for $\sim Q/\pi f$ and are weakly visible at the centre. Probing at
$r=0.9a$ with a $10^{-6}$ threshold and a minimum run time of 500 removed the plateau (the original results are kept in
`results/shortrun/`). An earlier hypothesis that a long run at the *same* probe would help was wrong; the diagnosis above
was established by the run-time study.

### Other Tier 1 checks

* **Mirror symmetry** (half domain) reproduces the full-domain run to $10^{-13}$ for TM and TE (integration tests).
* **Two independent total-cross-section estimates** (far-field integral and scattered flux) agree with each other and
  with the exact $\sigma_s$ at the levels in the table above; the flux normalisation `flux_norm_check` is 1.000.
* **Exact-solution tests** (unit tests): optical theorem, echo width integrating to the total width, interface-condition and boundary-solve checks, Rayleigh limit and exponent, PEC geometric-optics limit, truncation convergence, regression values.

## Tier 2 — periodic array workflow

### Verification against the exact slab (**COMPLETE**)

A homogeneous slab ($n=2$, $d=1$) has the closed-form Airy solution (`periodic.slab_rt`). The same periodic-cell
machinery (Bloch boundary, oblique source phase, flux normalisation, PML) is used for the cylinder array.

* Normal incidence, RMS error of $R$: 1.7e-2 (res 20), 7.4e-3 (30), 4.1e-3 (40), 1.8e-3 (60), 1.0e-3 (80);
  fitted order **2.03** (identical for TM and TE, as expected for a homogeneous layer).
* Oblique incidence at fixed $k_y=0.5, 1, 1.5$ (res 40), TM and TE: maximum $|R_\text{FDTD}-R_\text{exact}|\le1.1\times10^{-2}$.
* Energy residual $|1-R-T|\le1.2\times10^{-4}$ at normal incidence. At $k_y=1.5$ it reaches 1.2–1.6 % in the lowest
  frequencies of the scan ($f=0.30$–0.34, where $\sin\theta=k_y/k\approx0.8$, i.e. near-grazing incidence). This is a limitation of the
  finite-size source/monitor setup near grazing incidence; no conclusions in this repository rely on $\sin\theta>0.7$.

### Convergence of the resonance metrics (**COMPLETE**)

Lowest resonance of two representative arrays (`experiments/tier2/convergence_periodic.py`; 321-point windows):

| Case | Variation | $f_0$ range | $Q$ range |
|---|---|---|---|
| TM, Λ=3 | resolution 20→60 | 0.25470 → 0.25483 (5×10⁻⁴ relative) | 64.9 → 65.0 |
| TM, Λ=3 | margin 4–9, PML 3–6, monitor gap 2–5 | 0.25481 (identical to 5 digits) | 65.0 |
| TE, Λ=4 | resolution 20→60 | 0.19863 → 0.19870 (3.5×10⁻⁴) | 12.1 |
| TE, Λ=4 | margin, PML, monitor gap | 0.19869 | 12.1 |

The resonance frequency converges with resolution at the sub-10⁻³ level already at resolution 20; boundary and
monitor placement are irrelevant at the quoted precision. Resolution 40 is used for the studies.

### Independent Q estimates (**COMPLETE for the cases listed**)

Q from a coupled-mode line-shape fit of $R(f)$ is compared with the time-domain Harminv ringdown. TM, Λ=3: 59.7 (fit,
wide window) / 65.0 (fit, fine window) vs. **62.7** (Harminv): agreement within 7 %. Broad resonances ($Q\lesssim10$)
are not well defined by any single-Lorentzian model (Fano background): the line-shape fit and the FWHM estimate differ by up
to a factor of ~2 there, and Harminv does not detect such low-Q modes at all. **Q values below ~10 are therefore indicative only.**

### Known limitations of Tier 2 verification

* No second independent solver (e.g. RCWA) is included; verification rests on the exact slab, energy conservation,
  convergence and the ringdown cross-check (**PLANNED**: RCWA).
* The coarse 551-point survey under-samples lines narrower than ≈ 4×10⁻⁴ in frequency (Q ≳ 700) and some short-pitch TM runs did not
  decay within the time limit; those points are flagged in `docs/results.md` and excluded from conclusions.
* Lossless dielectric only; two-dimensional (z-invariant) geometry only.

## How to re-run

See the README. Fast suite without MEEP: `pytest -m "not meep"`; MEEP integration tests: `pytest -m meep`; stored-result checks:
`pytest -m results`.
