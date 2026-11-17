# Surface glint versus off-nadir angle and wind speed

> **Date:** 2026-11-16
> **Purpose:** surface glint (bin 0) as a function of the off-nadir viewing angle and of the
> wind speed, for the glint panel of the poster.
> **Code:** `v1.0-aisam2026` (the `.DAT` headers report `code = v0.8-receiver-dt`: same source
> code, tagged `v1.0-aisam2026` afterwards).
> **Duration:** 30 runs, 8 in parallel (4 physical cores, 8 threads): 68 s in total,
> 12–22 s per run.

## Command
From the repository root, after `make`:
```bash
runs/2026-11-16_glint_sweep/run_sweep.sh runs/2026-11-16_glint_sweep 1000000 8
```
The script creates `aXX_UYY.Y/input.nml` from `input/reference.nml` and runs
`mc_water input.nml > run.log 2>&1` in each subdirectory, 8 at a time (`xargs -P 8`).
Each subdirectory also contains `duration.txt` (seconds).

Figures, from the repository root:
```bash
python3 plot/fig_glint.py runs/2026-11-16_glint_sweep
```
The script writes `fig_glint_u10.pdf` and `fig_glint_alpha.pdf`; a copy of both is in this
directory.

## Input
`input.nml` = `input/reference.nml` (FOV = beam_div = 30 µrad, $H$ = 450 km,
$r_{tel}$ = 0.5 m, $c$ = 0.15 m$^{-1}$, $\omega$ = 0.85, $g$ = 0.92, seed 1,
`flat_sea = .false.`), with:
- `theta_v` ∈ {0, 6, 12, 20, 30} deg;
- `u10` ∈ {0, 2, 5, 7.5, 10, 15} m/s, giving $\sigma^2 = 0.003 + 0.00512\,U_{10}$ =
  0.003, 0.01324, 0.0286, 0.0414, 0.0542, 0.0798. At $U_{10} = 0$ the sea is slightly
  rough, not flat;
- `nphotons` = $10^6$.

## Results: glint (bin 0 of `SIGNAL_PROFILE_Z.DAT`)
| $\alpha$ \ $U_{10}$ [m/s] | 0 | 2 | 5 | 7.5 | 10 | 15 |
|---|---|---|---|---|---|---|
| 0° | 1.785e-12 | 4.044e-13 | 1.873e-13 | 1.293e-13 | 9.871e-14 | 6.708e-14 |
| 6° | 4.564e-14 | 1.784e-13 | 1.293e-13 | 1.007e-13 | 8.186e-14 | 5.939e-14 |
| 12° | 5.492e-19 | 1.424e-14 | 4.124e-14 | 4.641e-14 | 4.590e-14 | 4.070e-14 |
| 20° | 1.429e-31 | 2.199e-17 | 2.196e-15 | 6.351e-15 | 1.033e-14 | 1.537e-14 |
| 30° | 1.527e-60 | 7.249e-24 | 2.499e-18 | 6.343e-17 | 3.243e-16 | 1.583e-15 |

- **No zero glint.** The glint is computed analytically from the Cox–Munk slope
  distribution, so it is always positive. In the figures the y axis starts at $10^{-26}$:
  the two smallest values (20° and 30° at 0 m/s) are drawn on the axis as downward triangles.
- **At nadir (0°)** the glint decreases monotonically with wind speed, from 1.79e-12 at 0 m/s
  to 6.7e-14 at 15 m/s, i.e. a factor of 27.
- **At 6°** the glint peaks at 2 m/s and then decreases.
- **At 12°** it grows from 0 to 7.5 m/s (5.5e-19 → 4.64e-14), has a flat maximum between
  7.5 and 10 m/s and decreases by 12% at 15 m/s.
- **At 20° and 30°** it grows monotonically up to 15 m/s.
- **Glint versus subsurface signal at 12° and 7.5 m/s:** bin 0 (4.641e-14) is **45.7** times
  the maximum of the subsurface bins (1.016e-15, bin 1).

Figures: `fig_glint_u10.pdf` (glint versus $U_{10}$, one curve per angle) and
`fig_glint_alpha.pdf` (glint versus $\alpha$, one curve per wind speed).
