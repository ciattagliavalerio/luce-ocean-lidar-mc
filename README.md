# LUCE ocean lidar Monte Carlo

Semi-analytic Monte Carlo model of the coupled air–sea system for a spaceborne ocean lidar,
developed for the geometry of the LUCE mission (ASI–NASA): 450 km altitude, 12° off nadir, 532 nm.
The model simulates the time-based profile that a LUCE-like lidar records from the sea, the
corresponding profile in real depth, and the surface glint. The slant-view geometry, the Gaussian
beam, the wind-roughened sea surface and multiple scattering in water are taken into account.

This repository contains the code, the input files, the runs and the scripts used to produce every
figure of the poster presented at the AISAM Students' Conference 2026 (Festival della Meteorologia,
Rovereto, 21 November 2026). Version: `v1.0-aisam2026`.

**Author:** Valerio Ciattaglia (Physics, Università degli Studi di Roma Tor Vergata)

<!-- TODO: poster preview image, e.g. ![Poster](poster/poster_preview.png) -->

## The model

Each photon goes through six steps:

1. The photon is launched from the satellite with a Gaussian angular distribution of the beam
   (standard deviation `beam_div/4` per angular component).
2. At the air–sea interface the photon is refracted through a wind-roughened surface, whose facet
   slopes follow the Cox–Munk distribution, with Smith shadowing.
3. In water, free paths are sampled from the beam attenuation coefficient $c$, and scattering is
   described by the single-scattering albedo $\omega$ and the Henyey–Greenstein phase function
   (asymmetry $g$).
4. Photons outside the field of view (FOV) of the receiver keep travelling, but their
   contributions are recorded only while they are inside the FOV.
5. At every collision the probability that the photon reaches the telescope is computed
   analytically and scored (local estimate). The surface glint is computed with the same
   semi-analytic approach.
6. The contributions are binned both by real depth $z$ and by time-based depth $z_t$, separately
   for each scattering order. The time-based depth is the time of flight converted to depth and
   projected on the vertical with $\cos\alpha_r$ (Shangguan et al. 2025, Eq. 24 and 27).

## Repository structure

```
├── src/            Fortran 90 source: MC.f90 (main program), mc_module.f90 (procedures)
├── Makefile        build of the executable mc_water
├── input/          input templates (Fortran namelists)
├── runs/           one directory per run: input, output, log and README
├── plot/           Python scripts of the figures and common style
├── poster/         poster source (LaTeX, tikzposter) and PDF
├── LICENSE         MIT (code)
└── CITATION.cff    how to cite this repository
```

## Requirements

- `gfortran` and `make`;
- Python 3 with `numpy` and `matplotlib`, plus `qrcode` for the QR code (`pip install -r plot/requirements.txt`);
  the figures use the Lato font (if missing, matplotlib falls back to its default font);
- to compile the poster: XeLaTeX, `latexmk`, the packages `tikzposter`, `tcolorbox`, `qrcode`,
  `fontspec`, `mathastext`, and the Lato font.

## Quick start

```bash
make                                   # builds ./mc_water
mkdir -p runs/my_run && cp input/reference.nml runs/my_run/input.nml
cd runs/my_run
../../mc_water input.nml > run.log 2>&1
```

Every variable of the input file is required and documented in the comments of
`input/reference.nml`. Running time: about 1 min per $10^6$ photons on one core.

Output files (in the run directory):

| File | Content |
|---|---|
| `SIGNAL_PROFILE_Z.DAT` | signal per metre versus real depth $z$; bin 0 = surface (glint) |
| `SIGNAL_PROFILE_S.DAT` | signal per metre versus time-based depth $z_t$ |
| `SCATTERING_COUNTS.DAT` | collision statistics by scattering order |
| `input_used.nml` | the input actually used, with the random seed |

Each `.DAT` file starts with a header (`# key = value`) reporting the code version, the input
parameters and the output grid; then one row per bin: bin index, bin-centre depth [m], total
signal and signal of scattering orders 1–5.

## Reproducing the poster

From the repository root, after `make`:

| Poster figure | Run | Run command (from the run directory) | Figure command |
|---|---|---|---|
| Validation against ICESat-2 (`fig_icesat.pdf`) | `runs/2026-11-17_icesat2_c0138` | `../../mc_water input.nml > run.log 2>&1` | `python3 plot/fig_icesat.py` |
| Profiles by scattering order (`fig_orders_z.pdf`, `fig_orders_t.pdf`) | `runs/2026-11-17_luce_orders` | `../../mc_water input.nml > run.log 2>&1` | `python3 plot/fig_orders.py` |
| Glint versus wind speed (`fig_glint_u10.pdf`) | `runs/2026-11-16_glint_sweep` | `runs/2026-11-16_glint_sweep/run_sweep.sh runs/2026-11-16_glint_sweep 1000000 8` (from the root) | `python3 plot/fig_glint.py` |
| QR code of this repository (`qr_repo.png`) | — | — | `python3 plot/qr_repo.py` |

The figures are written to `poster/figures/`. Then:

```bash
python3 plot/style.py      # poster/palette.tex (only if the style changes)
cd poster && latexmk poster_LUCE.tex
```

All runs use seed 1, so they are reproducible bit by bit with the same compiler. Each run
directory has a `README.md` with purpose, parameters, duration and results.

**Note on the ICESat-2 figure.** The comparison uses curves digitized from Shangguan et al.
(2025), Fig. 4b, which are not distributed with this repository. Without them,
`plot/fig_icesat.py` draws only the profiles of this work; to add the published curves, see the
instructions at the top of the script.

## References

1. Cox, C., & Munk, W. (1954). *J. Opt. Soc. Am.* 44, 838.
2. Nakajima, T., & Tanaka, M. (1983). *J. Quant. Spectrosc. Radiat. Transfer* 29, 521.
3. Leathers, R. A., et al. (2004). *Monte Carlo Radiative Transfer Simulations for Ocean Optics:
   A Practical Guide*. NRL/MR/5660-04-8819.
4. Shangguan, M., et al. (2025). *Remote Sens. Environ.* 329, 114910.
5. Rozanov, V. V., et al. (2014). *J. Quant. Spectrosc. Radiat. Transfer* 133, 13.
6. Smith, B. G. (1967). *IEEE Trans. Antennas Propag.* 15, 668.
7. Henyey, L. G., & Greenstein, J. L. (1941). *Astrophys. J.* 93, 70.

<!-- TODO: add DOIs. -->

## License

- Code (`src/`, `plot/`, `Makefile`, scripts): [MIT](LICENSE).
- Poster (`poster/`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
  The university logo in `poster/figures/logos/` belongs to its owner and is not covered by
  these licenses.

## Citation

See [`CITATION.cff`](CITATION.cff), or use the "Cite this repository" button on GitHub.
