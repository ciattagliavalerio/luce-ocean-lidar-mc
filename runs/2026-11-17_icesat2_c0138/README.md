# ICESat-2 case of Shangguan et al. (2025), Fig. 4b

> **Date:** 2026-11-17
> **Purpose:** validation of the model against the ICESat-2 case of Shangguan et al. (2025),
> Fig. 4b; source of `fig_icesat.pdf` on the poster.
> **Code:** `v1.0-aisam2026`.
> **Duration:** 235 s on one core.

## Command
From this directory, after `make` in the repository root:
```bash
../../mc_water input.nml > run.log 2>&1
```

## Input
`input.nml` (template: `input/icesat2_shangguan2025.nml`):
ICESat-2 geometry, $\alpha$ = 2°, FOV = 88.4 µrad, beam_div = 35 µrad, $H$ = 500 km,
$r_{tel}$ = 0.4 m; flat sea without internal reflection at the surface (as in Shangguan et al.
2025, Sec. 4.1); $c$ = 0.1376 m⁻¹, $\omega$ = 0.5569, $g$ = 0.919; $10^8$ photons, seed 1.

Output grid from the receiver sampling interval (`&receiver`):
$$\cos\alpha_r = \sqrt{1 - (\sin\theta_v/n_w)^2} = 0.999655666,\qquad
\Delta t = \frac{2\,n_w\,dz}{c\,\cos\alpha_r} = 8.875861\ \mathrm{ns}$$
with $\theta_v = 2^\circ$, $n_w = 1.33$, $dz = 1$ m and $c = 299792458$ m/s. With
`dt_bin` = 8.876 ns, $dz$ = 1.0000156 m; with `n_bins` = 200, $z_{max}$ = 200.003 m.

## Comparison with Shangguan et al. (2025), Fig. 4b
Reference values were digitized from Shangguan et al. (2025), Fig. 4b; the digitized curves are
not distributed with this repository. Profiles are normalized to the larger of the two maxima.
"60 dB" = fit from the surface to the depth where the profile drops to $10^{-6}$ ($z_d$ for the
time-based profile, $z^{gt}$ for the real-depth profile).

| Quantity | This run | Shangguan Fig. 4b | Difference |
|---|---|---|---|
| $K_t$, 60 dB [m⁻¹] | 0.0722 ± 0.0005 | 0.0719 | +0.4% |
| $K_z$, 60 dB [m⁻¹] | 0.1049 ± 0.0002 | 0.1051 | −0.2% |
| $K_z/K_t$, 60 dB | 1.45 | 1.46 | −0.6% |
| $z_d$ [m] | 93.5 | 94.8 | −1.4% |
| $z^{gt}$ [m] | 67.5 | 66.7 | +1.2% |
| $K_t$, 5–60 m [m⁻¹] | 0.0675 ± 0.0005 | 0.0678 | −0.5% |
| $K_z$, 5–60 m [m⁻¹] | 0.1048 ± 0.0002 | 0.1048 | 0.0% |
| $K_z/K_t$, 5–60 m | 1.55 | 1.55 | 0% |
| Surface bin | 0 | — | — |

**Point-by-point deviation from the digitized curves** (down to 100 m;
$\Delta = \log_{10}$ ours $-$ $\log_{10}$ digitized, our profile interpolated at the digitized depths):

| Profile | Median of $\lvert\Delta\rvert$ | Mean of $\Delta$ | Maximum of $\lvert\Delta\rvert$ |
|---|---|---|---|
| Real depth | 0.052 (factor 1.13) | +0.046 | 0.114 (factor 1.30) |
| Time-based | 0.073 (factor 1.18) | +0.034 | 0.552 (factor 3.6) |

The attenuation coefficients agree within 1%. The positive mean deviation means that, with this
normalization, our profiles lie slightly above the published ones. The largest deviations of the
time-based profile are all deeper than 89 m (the largest at 96–97 m), in the tail near $10^{-6}$, where
both the Monte Carlo profile and the digitization are least accurate.

The surface bin is 0: on a flat sea at $\alpha$ = 2° the specular reflection does not reach the
telescope, so there is no glint.

## Figure
`fig_icesat.pdf` on the poster (script `plot/fig_icesat.py`).
