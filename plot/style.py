"""
style.py
--------
Single source of style for the LUCE poster: palette, fonts, sizes.

Every figure script imports this module and never defines colours or
sizes on its own. Running it as a program writes poster/palette.tex, so
that the LaTeX source uses exactly the same colours:

    python3 plot/style.py

The font is Lato: if it is not installed, matplotlib falls back to its
default font and the figures look slightly different.

Figures are produced at their final printed size (scale 1:1 on the A0
sheet), so the font sizes below are the real ones on paper.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. Palette (Okabe-Ito, colour-blind safe). One colour, one quantity.
# ----------------------------------------------------------------------
COLORS = {
    "real":      "#0072B2",   # Real-depth profile, K_z.
    "time":      "#D55E00",   # Time-based profile, K_t.
    "order1":    "#009E73",   # Single scattering.
    "order2":    "#E69F00",   # Second order.
    "order3":    "#56B4E9",   # Third order.
    "order4":    "#CC79A7",   # Orders >= 4.
    "literature": "#000000",  # Published data (markers, no line).
    "limit":     "#666666",   # Analytic limits (dashed).
    "frame":     "#12355B",   # Poster frames, tabs, title.
    "keyfill":   "#FFF4CC",   # Key-message boxes.
    "panelfill": "#F3F6FA",   # Inner panels.
    "grid":      "#D9D9D9",   # Plot grid.
    "text":      "#1A1A1A",   # Body text.
    "placeholder": "#9A9A9A", # Placeholders.
}
# Ordered series (angles, wind speeds): the same Okabe-Ito colours as the
# profile graphs, in order blue, sky blue, green, orange, vermilion.
SERIES_COLORS = [COLORS["real"], COLORS["order3"], COLORS["order1"],
                 COLORS["order2"], COLORS["time"], COLORS["order4"]]

# Colours of the steps of the photon cycle (Methods cells and big drawing).
# A separate group: green, electric blue, indigo, magenta, crimson, so that a step is
# never mistaken for a quantity of the data palette above. Vivid colours,
# chosen so that every pair (steps among themselves and steps against the
# data colours) stays apart also in simulated deuteranopia and protanopia
# (minimum CIELAB distance about 12), and white bold titles on them keep a
# contrast ratio of at least 3.2. Step 6 (profiles) uses "real" and "time".
STEP_COLORS = {
    "step1": "#23A443",   # Launch: Gaussian beam (bright green).
    "step2": "#006AFF",   # Air-sea interface: waves, facet normal (electric blue).
    "step3": "#452797",   # In-water walk: path, collisions, theta (indigo).
    "step4": "#B61C84",   # FOV check: receiver cone, photon outside (magenta).
    "step5": "#B01F24",   # Local estimate: scores and glint (crimson).
}

# ----------------------------------------------------------------------
# 2. Sizes
# ----------------------------------------------------------------------
CM = 1.0 / 2.54               # Inches per centimetre.
FIG_W_CM = 18.6               # Standard graph width [cm].
FIG_H_CM = 13.5               # Standard graph height [cm].
FONT_PT = 24                  # Axis labels and ticks [pt].
LEGEND_PT = 28                # Legends, larger than the labels [pt].
LINE_PT = 3.0                 # Line width [pt].
ORDER_LINE_PT = 1.8           # Thinner lines for the scattering orders [pt].
HALO_PT = 2.5                 # White halo around the main curves, drawn on top [pt].
MARKER_PT = 10                # Marker size [pt].
FONT_FAMILY = "Lato"

# ----------------------------------------------------------------------
# 3. Labels and units: the same quantity is always written the same way
# ----------------------------------------------------------------------
LABELS = {
    "depth":   r"Depth $z$, $z_t$ [m]",
    "depth_z": r"Depth $z$ [m]",          # Real-depth profile only.
    "depth_t": r"Depth $z_t$ [m]",        # Time-based profile only.
    "signal":  "Normalized signal",
    "K":       r"$K$ [m$^{-1}$]",
    "fov":     r"FOV [$\mu$rad]",
    "order1":  "Order 1",
    "order2":  "Order 2",
    "order3":  "Order 3",
    "order4":  r"Orders $\geq$ 4",
    "u10":     r"$U_{10}$ [m s$^{-1}$]",
    "alpha":   r"$\alpha$ [deg]",
    "glint":   "Glint signal",
    "real":    "Real depth",
    "time":    "Time-based",
    "Kz":      r"$K_z$",
    "Kt":      r"$K_t$",
}


def apply():
    """Set the matplotlib defaults shared by all the figures."""
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": FONT_PT,
        "axes.labelsize": FONT_PT,
        "axes.titlesize": FONT_PT,
        "xtick.labelsize": FONT_PT,
        "ytick.labelsize": FONT_PT,
        "legend.fontsize": LEGEND_PT,
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_FAMILY,
        "mathtext.it": FONT_FAMILY + ":italic",
        "mathtext.bf": FONT_FAMILY + ":bold",
        "mathtext.cal": FONT_FAMILY,
        "mathtext.fallback": "cm",
        "lines.linewidth": LINE_PT,
        "lines.markersize": MARKER_PT,
        "axes.linewidth": 1.8,
        "axes.edgecolor": COLORS["frame"],
        "axes.labelcolor": COLORS["text"],
        "text.color": COLORS["text"],
        "xtick.color": COLORS["text"],
        "ytick.color": COLORS["text"],
        "xtick.major.width": 1.8,
        "ytick.major.width": 1.8,
        "xtick.major.size": 8,
        "ytick.major.size": 8,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "axes.grid": True,
        "grid.color": COLORS["grid"],
        "grid.linewidth": 1.0,
        "axes.axisbelow": True,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.borderaxespad": 0.3,
        "legend.labelspacing": 0.25,
        "pdf.fonttype": 42,
        "savefig.transparent": True,
    })


def halo():
    """White outline that keeps a main curve readable over the others."""
    from matplotlib import patheffects as pe
    return [pe.Stroke(linewidth=LINE_PT + 2 * HALO_PT, foreground="white"), pe.Normal()]


def new_figure():
    """Figure and axes of the standard size."""
    apply()
    fig, ax = plt.subplots(figsize=(FIG_W_CM * CM, FIG_H_CM * CM),
                           constrained_layout=True)
    return fig, ax


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSTER = os.path.join(ROOT, "poster")         # Poster source directory.
FIGURES = os.path.join(POSTER, "figures")     # Where the figures are written.


def save(fig, name):
    """Write the figure as vector PDF in poster/figures/."""
    os.makedirs(FIGURES, exist_ok=True)
    out = os.path.join(FIGURES, name + ".pdf")
    fig.savefig(out)
    plt.close(fig)
    print("written", out)


def write_palette_tex():
    """Write poster/palette.tex with the same colours for the LaTeX source."""
    out = os.path.join(POSTER, "palette.tex")
    with open(out, "w") as f:
        f.write("% Generated by plot/style.py: do not edit by hand.\n")
        for name, hexa in list(COLORS.items()) + list(STEP_COLORS.items()):
            f.write("\\definecolor{luce%s}{HTML}{%s}\n"
                    % (name, hexa.lstrip("#")))
        f.write("\\newcommand{\\figwcm}{%.1f}\n" % FIG_W_CM)
        f.write("\\newcommand{\\fighcm}{%.1f}\n" % FIG_H_CM)
    print("written", out)


if __name__ == "__main__":
    write_palette_tex()
