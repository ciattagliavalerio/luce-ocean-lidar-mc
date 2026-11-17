"""
fig_icesat.py - ICESat-2 case against Shangguan et al. (2025), Fig. 4b.

Data: runs/2026-11-17_icesat2_c0138 (this work, code v1.0-aisam2026).
Optional overlay: curves digitized from Shangguan et al. (2025), Fig. 4b,
which are NOT distributed with this repository. To draw them, put two
CSV files in data/shangguan2025_fig4b/ (or in the directory given by the
environment variable SHANGGUAN_FIG4B):
    fig4b_real.csv  (real-depth profile)
    fig4b_time.csv  (time-based profile)
each with two comma-separated columns: depth [m], log10 of the normalized
signal. Without them only the profiles of this work are drawn.
Usage, from the repository root: python3 plot/fig_icesat.py
Writes: poster/figures/fig_icesat.pdf
"""

import os
import numpy as np
from run_data import ROOT, RUNS, load_profiles
import style as st

RUN = os.path.join(RUNS, "2026-11-17_icesat2_c0138")
DIGIT = os.environ.get("SHANGGUAN_FIG4B",
                       os.path.join(ROOT, "data", "shangguan2025_fig4b"))


def thin(depth, step):
    """Indices of about one digitized point every `step` metres."""
    idx, last = [], -1e9
    for i in np.argsort(depth):
        if depth[i] - last >= step:
            idx.append(i)
            last = depth[i]
    return np.array(idx)


def main():
    par, dz, sz, dt, s_t, _ = load_profiles(RUN)
    norm = max(sz[:, 0].max(), s_t[:, 0].max())
    files = [os.path.join(DIGIT, "fig4b_%s.csv" % k) for k in ("real", "time")]
    have_digit = all(os.path.isfile(f) for f in files)
    if not have_digit:
        print("WARNING: digitized curves of Shangguan et al. (2025) not found in",
              DIGIT, "- drawing only the profiles of this work.")
    fig, ax = st.new_figure()
    mz, mt = sz[:, 0] > 0, s_t[:, 0] > 0
    ax.semilogy(dz[mz], sz[mz, 0] / norm, color=st.COLORS["real"],
                label=st.LABELS["real"], zorder=5, path_effects=st.halo())
    ax.semilogy(dt[mt], s_t[mt, 0] / norm, color=st.COLORS["time"],
                label=st.LABELS["time"], zorder=5, path_effects=st.halo())
    # Published points: open and sparse, on top of the model curves.
    kw = dict(color=st.COLORS["literature"], ls="none", ms=9, mfc="none",
              mew=1.6, zorder=6)
    for k, f in enumerate(files if have_digit else []):
        pts = np.loadtxt(f, delimiter=",")
        sel = thin(pts[:, 0], 5.0)
        ax.semilogy(pts[sel, 0], 10.0 ** pts[sel, 1], marker="o", **kw,
                    label="Shangguan et al." if k == 0 else None)
    ax.set_xlim(0, 100)
    ax.set_ylim(1e-15, 3)                     # Room for the legend.
    ax.set_yticks([1e-12, 1e-9, 1e-6, 1e-3, 1])
    ax.set_xlabel(st.LABELS["depth"])
    ax.set_ylabel(st.LABELS["signal"])
    ax.legend(loc="lower left")
    st.save(fig, "fig_icesat")


if __name__ == "__main__":
    main()
