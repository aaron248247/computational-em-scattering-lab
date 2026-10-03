# Periodic dielectric-cylinder array

The periodic calculation uses the same dielectric cylinder as the single-cylinder scattering problem, but places one cylinder in each unit cell of an infinite one-dimensional lattice. Varying the lattice pitch introduces coupling between neighboring cylinders and changes the resonance spectrum.

## 1. Geometry, polarization, and observables

Cylinders of radius $a$ and permittivity $\varepsilon_r$ are centered at

$$
y_j=j\Lambda,
$$

with their axes along $z$ and centers lying on $x=0$. The incident plane wave travels primarily along $+x$ and may have transverse wavenumber $k_y$. Normal incidence corresponds to $k_y=0$, and

$$
\sin\theta=\frac{k_y}{k},\qquad k=\frac{2\pi f}{c}.
$$

Polarization follows the same convention as the single-cylinder calculation:

- TM: $E_z$;
- TE: $H_z$.

The main observables are reflectance $R$ and transmittance $T$, obtained from Poynting flux through planes on either side of the array and normalized by an empty-cell reference run. For a lossless structure,

$$
A=1-R-T
$$

should vanish apart from numerical error.

## 2. Bloch periodicity and diffraction orders

For a structure periodic in $y$, Bloch's theorem gives

$$
F(x,y)=e^{ik_y y}u(x,y),
$$

where $u(x,y)$ is periodic with period $\Lambda$. A single unit cell with Bloch-periodic boundaries can therefore represent the infinite array.

Outside the array, the field can be written as a sum of diffraction orders with transverse wavenumbers

$$
k_{y,m}=k_y+\frac{2\pi m}{\Lambda}.
$$

Order $m$ propagates when

$$
|k_{y,m}|<k
$$

and is evanescent otherwise.

At normal incidence, only the zeroth order propagates when $f<c/\Lambda$. New propagating orders appear at the Rayleigh-Wood anomaly frequencies

$$
f_m=\frac{|k_y+2\pi m/\Lambda|}{2\pi}.
$$

These frequencies are computed by `periodic.wood_anomaly_frequencies`. The distinction between Rayleigh-type order openings and resonance-type Wood anomalies is discussed by [Hessel and Oliner (1965), *A New Theory of Wood's Anomalies on Optical Gratings*](https://doi.org/10.1364/AO.4.001275).

## 3. Resonant response

The reflectance spectrum contains contributions from both local and collective resonances.

Local, Mie-like resonances originate from the same cylinder modes present in the isolated-scatterer problem. Their frequencies and linewidths change when neighboring cylinders are brought close enough to interact.

The array also supports collective Bloch modes. Modes above the light line can couple to radiating diffraction orders and therefore have finite radiative linewidth; guided-resonance behavior of this type is analyzed by [Fan and Joannopoulos (2002), *Analysis of Guided Resonances in Photonic Crystal Slabs*](https://doi.org/10.1103/PhysRevB.65.235112). Modes below the light line are guided and do not radiate into the surrounding medium.

Near an isolated spectral feature, the reflected amplitude is modeled with the temporal coupled-mode form developed for Fano resonances by [Fan, Suh, and Joannopoulos (2003), *Temporal Coupled-Mode Theory for the Fano Resonance in Optical Resonators*](https://doi.org/10.1364/JOSAA.20.000569)

$$
r(f)=r_b+\frac{A}{f-f_0+i\gamma},\qquad Q=\frac{f_0}{2\gamma},
$$

where $r_b$ is a slowly varying background, $f_0$ is the resonance frequency, and $\gamma$ is the radiative half-width. Interference between the background and resonant terms produces an asymmetric Fano line shape. For a lossless mirror-symmetric layer, a strongly coupled resonance can approach $R=1$.

The fitted quantities are computed in `em_rf_lab.resonance` and summarized by `em_rf_lab.tier2_analysis`.

## 4. Symmetry and bound states in the continuum

For a centered cylinder, the unit cell is mirror-symmetric under $y\to-y$ and $x\to-x$. At $k_y=0$, a mode whose parity is incompatible with that of the radiating zeroth-order plane wave cannot couple to that channel. Such a mode can remain localized even though its frequency lies in the radiation continuum, producing a symmetry-protected bound state in the continuum (BIC); a broad review of this mechanism is given by [Hsu et al. (2016), *Bound States in the Continuum*](https://doi.org/10.1038/natrevmats.2016.48).

Moving away from $k_y=0$ breaks the symmetry and introduces radiative loss. Perturbation theory commonly gives

$$
\frac{1}{Q}\propto k_y^2
$$

near the $\Gamma$ point for a symmetry-protected BIC.

The present calculations do not establish this scaling. The Bloch-wavenumber scan is preliminary, and a dedicated near-$\Gamma$ study with mode-symmetry analysis would be required before making a BIC claim.

The large $Q$ values seen when varying pitch at normal incidence are instead consistent with an accidental near-decoupling of a mode from the zeroth radiation channel. That interpretation is treated as a hypothesis rather than a demonstrated mechanism; see [`05_results.md`](05_results.md).

## 5. Checks used for the periodic calculation

The periodic-cell implementation is checked in several ways:

- A homogeneous dielectric slab is compared with the exact Airy/transfer-matrix solution implemented in `periodic.slab_rt`.
- Energy conservation is checked through $R+T=1$ for lossless structures.
- Rayleigh-Wood anomaly positions are calculated directly from the grating equation.
- Spectral $Q$ from line-shape fitting is compared with a time-domain Harminv ringdown estimate where the resonance is sufficiently narrow and long-lived.
- Mirror symmetry and reciprocity are checked in the computed spectra.
- In the large-pitch limit, the lowest collective resonance should approach the corresponding isolated-cylinder resonance. That limit has not yet been completed.
