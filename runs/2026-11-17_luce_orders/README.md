# LUCE case: scattering orders

> **Date:** 2026-11-17
> **Purpose:** real-depth and time-based profiles of the LUCE-like configuration, split by
> scattering order; source of `fig_orders_z.pdf` and `fig_orders_t.pdf` on the poster.
> **Code:** `v1.0-aisam2026`.
> **Duration:** 554 s on one core (run in parallel with `2026-11-17_icesat2_c0138`).

## Command
From this directory, after `make` in the repository root:
```bash
../../mc_water input.nml > run.log 2>&1
```

## Input
`input.nml` = `input/reference.nml` with `nphotons` = $10^8$ (instead of $10^7$):
$\alpha$ = 12°, FOV = beam_div = 30 µrad, $H$ = 450 km, $r_{tel}$ = 0.5 m,
$U_{10}$ = 7.5 m/s (rough sea), internal reflection at the surface on, $c$ = 0.15 m⁻¹,
$\omega$ = 0.85, $g$ = 0.92, seed 1.

Output grid from the receiver sampling interval (`&receiver`):
$$\cos\alpha_r = \sqrt{1 - (\sin\theta_v/n_w)^2} = 0.987705738,\qquad
\Delta t = \frac{2\,n_w\,dz}{c\,\cos\alpha_r} = 8.983247\ \mathrm{ns}$$
with $\theta_v = 12^\circ$, $n_w = 1.33$, $dz = 1$ m and $c = 299792458$ m/s. With
`dt_bin` = 8.983 ns, $dz$ = 0.9999725 m; with `n_bins` = 100, $z_{max}$ = 99.997 m.

## Results
Slope-method fit between 5 and 60 m, $E \propto e^{-2Kz}$, on the subsurface bins.

| Quantity | Value |
|---|---|
| $K_z$ [m⁻¹] | 0.1076 ± 0.0003 |
| $K_t$ [m⁻¹] | 0.0519 ± 0.0004 |
| $K_z/K_t$ | 2.07 |
| Crossing depth of real and time-based profiles | 7.5 m |
| First bin, share of orders 1…5, real depth | 0.326, 0.171, 0.133, 0.106, 0.083 |
| First bin, share of orders 1…5, time-based | 0.869, 0.120, 0.010, 0.001, 0.000 |
| Surface bin (glint) | 4.640e-14 |
| Subsurface maximum | 1.028e-15 |
| Glint / subsurface maximum | 45.1 |

## Figures
`fig_orders_z.pdf` and `fig_orders_t.pdf` on the poster (script `plot/fig_orders.py`).
