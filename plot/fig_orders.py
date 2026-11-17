"""
fig_orders.py - LUCE-like profile split by scattering order.

Data: runs/2026-11-17_luce_orders (1e8 photons, code v1.0-aisam2026).
Subsurface bins only; both profiles are normalized by the larger of the
two maxima. Orders 1, 2, 3 are read from the run; orders >= 4 are the
total minus orders 1-3 (orders() of run_data.py).
Usage, from the repository root: python3 plot/fig_orders.py
Writes: poster/figures/fig_orders_z.pdf (real depth),
poster/figures/fig_orders_t.pdf (time-based depth).
"""

import os
import numpy as np
from run_data import RUNS, load_profiles, orders
import style as st

RUN = os.path.join(RUNS, "2026-11-17_luce_orders")
KEYS = ("order1", "order2", "order3", "order4")


def draw(depth, sig, norm, total_key, name):
    """Total in the colour of its profile, thick with a white halo; orders
    dashed, thinner and on top, so that both stay visible where they overlap.

    Only the total enters the legend: the legend of the orders is shared
    by the two graphs and drawn under them in the poster (\\orderlegend).
    """
    tot, *by_order = orders(sig)
    fig, ax = st.new_figure()
    m = tot > 0
    ax.semilogy(depth[m], tot[m] / norm, color=st.COLORS[total_key],
                label=st.LABELS[total_key], zorder=5, path_effects=st.halo())
    for o, key in zip(by_order, KEYS):
        m = o > 0
        ax.semilogy(depth[m], o[m] / norm, color=st.COLORS[key], ls="--",
                    lw=st.ORDER_LINE_PT, zorder=6)
    ax.set_xlim(0, 100)
    ax.set_ylim(1e-15, 3)
    ax.set_yticks([1e-12, 1e-9, 1e-6, 1e-3, 1])
    # Each graph has its own depth: z (real) or z_t (time-based).
    ax.set_xlabel(st.LABELS["depth_z" if total_key == "real" else "depth_t"])
    ax.set_ylabel(st.LABELS["signal"])
    ax.legend(loc="lower left")
    st.save(fig, name)


def main():
    par, dz, sz, dt, s_t, _ = load_profiles(RUN)
    norm = max(sz[:, 0].max(), s_t[:, 0].max())
    draw(dz, sz, norm, "real", "fig_orders_z")
    draw(dt, s_t, norm, "time", "fig_orders_t")
    # Shares of the orders in the first subsurface bin (0-1 m).
    for label, sig in (("real", sz), ("time", s_t)):
        tot, *by_order = orders(sig)
        print(label, "first bin, shares of orders 1, 2, 3, >=4:",
              np.round(np.array(by_order)[:, 0] / tot[0], 3))


if __name__ == "__main__":
    main()
