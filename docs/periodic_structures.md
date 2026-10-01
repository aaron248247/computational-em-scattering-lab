# Periodic structures: theory for the primary investigation (Tier 2)

The Tier 2 study applies the verified FDTD workflow to a **one-dimensional periodic array of dielectric
cylinders** (a minimal frequency-selective surface / dielectric metasurface). The scatterer is exactly the
cylinder that is validated against Mie theory in Tier 1, so the investigation asks how *periodic arrangement*
converts isolated-cylinder resonances into collective resonances of a layer.

## 1. Geometry, polarisation, observables

Cylinders (radius $a$, permittivity $\varepsilon_r$) centred on $x=0$ at $y_j=j\Lambda$, invariant along $z$,
illuminated by a plane wave travelling along $+x$ with transverse wavenumber $k_y$ ($k_y=0$ is normal
incidence; $\sin\theta=k_y/k$, $k=2\pi f/c$). Polarisation names follow the benchmark: **TM** = $E_z$
(*s*-polarised), **TE** = $H_z$ (*p*-polarised).

Observables are the total reflectance and transmittance (Poynting flux through transverse planes divided by
the incident flux of an empty-cell reference run), summed over all propagating orders, and the absorption
residual $A=1-R-T$, which **must vanish** for a lossless layer and is used as a convergence/validity
diagnostic.

## 2. Bloch theorem and diffraction orders

A structure periodic in $y$ admits fields $F(x,y)=e^{ik_yy}\,u(x,y)$ with $u$ periodic, which is how a
single unit cell with Bloch-periodic boundaries represents the infinite array. Outside the layer the field is a
sum of orders with transverse wavenumber $k_{y,m}=k_y+2\pi m/\Lambda$; order $m$ **propagates** if
$|k_{y,m}|<k$ and is evanescent otherwise. At normal incidence only the zeroth order propagates for
$f<c/\Lambda$, so the layer is a two-port with a single reflection and transmission channel; new orders
switch on at the **Rayleigh–Wood anomalies** $f_m=|k_y+2\pi m/\Lambda|/2\pi$ (implemented in
`periodic.wood_anomaly_frequencies`), where the response has a characteristic cusp because power is redistributed.

## 3. Resonances: where the spectral features come from

Two classes of resonance shape $R(f)$:

1. **Mie-like (local) resonances** of each cylinder — the interior modes verified in Tier 1 — shifted and
   broadened by coupling to the neighbouring cylinders through the lattice.
2. **Collective (lattice / guided-mode) resonances**, where the field is a Bloch mode of the whole array.
   Above the light line these modes radiate into the zeroth order and are *leaky* (finite quality factor);
   below the light line they are truly guided (infinite $Q$, invisible in a plane-wave experiment).

Near an isolated resonance, **coupled-mode theory** gives the reflection amplitude

$$
r(f)=r_b+\frac{A}{f-f_0+i\gamma},\qquad Q=\frac{f_0}{2\gamma},
$$

with $r_b$ a slowly varying background (whose interference with the resonance produces asymmetric Fano line
shapes), $\gamma$ the radiative half-width and $R=|r|^2$. For a lossless, mirror-symmetric layer a resonance
that couples to the zeroth order can drive $R\to1$ (full reflection), the operating principle of a
resonant-reflection filter. Quality factor sets the bandwidth trade-off: the fractional 90 % reflection band
scales like $1/Q$. Quantities extracted in `em_rf_lab.resonance` and `em_rf_lab.tier2_analysis`.

## 4. Symmetry selection rules and bound states in the continuum

The layer is mirror-symmetric in $y\to-y$ (for a centred cylinder) and $x\to-x$. At $k_y=0$ a mode whose
field pattern has the wrong parity cannot couple to the symmetric zeroth-order plane wave, so it is
**decoupled from the radiation continuum** despite lying inside it: a *symmetry-protected bound state in the
continuum* (BIC), with $Q=\infty$. Away from $k_y=0$ the symmetry is broken and the mode acquires a radiative
width, and from perturbation theory $1/Q\propto k_y^{2}$ near the $\Gamma$ point. Such modes are
invisible in the normal-incidence spectrum, appear as extremely narrow features at small oblique angle, and
show up in the ringdown (Harminv) analysis with enormous $Q$. Whether this behaviour is realised, and its
scaling exponent, is **measured** in the study (`docs/results.md`), not assumed.

## 5. Independent reference and sanity checks

* **Exact slab (Airy / transfer-matrix) solution** for a homogeneous film of the same thickness and index:
  validates the periodic-cell machinery (Bloch boundary, oblique source phase, flux normalisation, PML) with
  no patterning (`periodic.slab_rt`; closed-form normal-incidence check and Brewster condition in the tests).
* **Energy conservation** $R+T=1$ for the lossless array.
* **Wood-anomaly positions** are exact (grating equation) and must coincide with the onset of new orders.
* **Two independent estimators of $Q$**: line-shape fit of $R(f)$ and time-domain ringdown (Harminv).
* **Mirror symmetry** and **reciprocity** of the computed spectra.
* **Limit check**: as $\Lambda\to\infty$ the lowest collective resonance should approach the isolated-cylinder
  Mie resonance computed with the Tier 1 reference.
