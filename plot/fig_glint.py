"""
fig_glint.py - Glint against wind speed and off-nadir angle.

Sweep: runs/2026-11-16_glint_sweep. Every run is a subfolder whose
SIGNAL_PROFILE_Z.DAT header gives alpha and U10 (any folder names).

Usage, from the repository root: python3 plot/fig_glint.py [<sweep_dir>]
(default: runs/2026-11-16_glint_sweep).
Writes: poster/figures/fig_glint_u10.pdf, poster/figures/fig_glint_alpha.pdf
"""

import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from run_data import RUNS, read_header
import style as st


def read_sweep(sweep_dir):
    """Dictionary (alpha, U10) -> glint (bin 0 of the real profile)."""
    data = {}
    for path in glob.glob(os.path.join(sweep_dir, "**", "SIGNAL_PROFILE_Z.DAT"),
                          recursive=True):
        par = read_header(path)
        data[(par["alpha"], par["U10"])] = np.loadtxt(path)[0, 2]
    return data


# Lower end of the y axis: smaller values (low wind, large angle, down to
# 1e-60) are drawn on the axis as downward triangles ("below the scale").
Y_FLOOR = 1e-26


def draw(data, x_index, name, xlabel, fmt, legend_title, legend_loc):
    """One curve per value of the other variable, data-palette colours."""
    others = sorted({k[1 - x_index] for k in data})
    fig, ax = st.new_figure()
    for i, other in enumerate(others):
        pts = np.array(sorted((k[x_index], v) for k, v in data.items()
                              if k[1 - x_index] == other and v > 0))
        x = pts[:, 0]
        y = np.maximum(pts[:, 1], Y_FLOOR)
        low = pts[:, 1] < Y_FLOOR
        color = st.SERIES_COLORS[i % len(st.SERIES_COLORS)]
        ax.semilogy(x, y, color=color, label=fmt % other, marker="o",
                    markevery=list(np.flatnonzero(~low)), clip_on=False)
        ax.plot(x[low], y[low], "v", color=color, clip_on=False)
    ax.set_ylim(bottom=Y_FLOOR)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(st.LABELS["glint"])
    # Legend in the empty corner left by the low end of the axis.
    ax.legend(title=legend_title, loc=legend_loc, ncol=2,
              handlelength=1.0, columnspacing=0.6, borderaxespad=0.3)
    st.save(fig, name)


def main():
    sweep = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        RUNS, "2026-11-16_glint_sweep")
    data = read_sweep(sweep)
    draw(data, 1, "fig_glint_u10", st.LABELS["u10"], "%g°",
         r"$\alpha$", "lower right")
    draw(data, 0, "fig_glint_alpha", st.LABELS["alpha"], "%g",
         r"$U_{10}$ [m s$^{-1}$]", "lower left")


if __name__ == "__main__":
    main()
