"""
run_data.py
-----------
Reading of the run files and slope-method fit, shared by the figure
scripts. No colours or sizes here: see style.py.

Output files of mc_water (SIGNAL_PROFILE_Z.DAT, SIGNAL_PROFILE_S.DAT):
a header of '# key = value' lines, then one row per bin with
col 1 = bin index, col 2 = bin-centre depth [m], cols 3-8 = signal per
metre (total, then scattering orders 1-5). Bin 0 is the surface (glint)
and is excluded from the subsurface profiles. In SIGNAL_PROFILE_S.DAT the
time-based depth is already projected on the vertical with cos(alpha_r)
(Shangguan et al. 2025, Eq. 27).
"""

import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "runs")             # Simulation runs.
N_WATER = 1.33                                # Refractive index of water.


def read_header(path):
    """Return the '# key = value' header lines as a dict.

    Numeric values become floats; other values of identifier-like keys
    (e.g. code, flat_sea = T/F) are kept as strings.
    """
    params = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            body = line.lstrip("#").strip()
            if "=" in body and not body.startswith("iz") and not body.startswith("col"):
                key, val = body.split("=", 1)
                key, val = key.strip(), val.strip()
                if not val:
                    continue
                try:
                    params[key] = float(val.split()[0])
                except ValueError:
                    if key.replace("_", "").isalnum():
                        params[key] = val
    return params


def load_profiles(run_dir):
    """Header, depth axes, signals (total and orders, bins >= 1) and surface bin.

    Returns par, depth_z, sig_z, depth_t, sig_t, surface_z.
    """
    fz = os.path.join(run_dir, "SIGNAL_PROFILE_Z.DAT")
    fs = os.path.join(run_dir, "SIGNAL_PROFILE_S.DAT")
    z = np.loadtxt(fz)
    s = np.loadtxt(fs)
    return read_header(fz), z[1:, 1], z[1:, 2:], s[1:, 1], s[1:, 2:], z[0, 2]


def orders(sig):
    """Total, orders 1, 2, 3 and orders >= 4 (total minus orders 1-3)."""
    if sig.shape[1] < 4:
        sys.exit("ERROR: the run has maxsca < 3; these plots need orders 1-3.")
    tot = sig[:, 0]
    ms1, ms2, ms3 = sig[:, 1], sig[:, 2], sig[:, 3]
    return tot, ms1, ms2, ms3, tot - ms1 - ms2 - ms3


def fit_k(depth, sig, zmin, zmax):
    """K from E ~ exp(-2 K z) between zmin and zmax, with the fit error."""
    m = (depth >= zmin) & (depth <= zmax) & (sig > 0)
    p, cov = np.polyfit(depth[m], np.log(sig[m]), 1, cov=True)
    return -p[0] / 2.0, np.sqrt(cov[0, 0]) / 2.0
