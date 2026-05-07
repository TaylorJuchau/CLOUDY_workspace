"""
Optional helper: estimate model AGN fraction at optical wavelengths.

Run after a CIGALE AGN fit. Point this script at a best-model SED FITS file,
for example:

    python estimate_optical_agn_fraction.py out/NGC1068_best_model.fits

CIGALE's exact component-column names vary by version and module. This script
prints the available columns and then tries to identify AGN and total model
luminosity columns automatically. If automatic matching fails, edit AGN_KEYS and
TOTAL_KEYS below.
"""

from __future__ import annotations

import sys
import numpy as np
from astropy.table import Table

AGN_KEYS = ["agn", "fritz", "skirtor", "torus"]
TOTAL_KEYS = ["total", "Fnu", "L_lambda_total", "luminosity"]
REST_WAVELENGTHS_MICRON = [0.44, 0.55, 0.80]


def find_column(cols, keys, avoid=()):
    for c in cols:
        lc = c.lower()
        if any(k.lower() in lc for k in keys) and not any(a.lower() in lc for a in avoid):
            return c
    return None


def interp_at(wave, y, lam_micron):
    return np.interp(lam_micron, wave, y, left=np.nan, right=np.nan)


def main(path):
    tab = Table.read(path)
    cols = list(tab.colnames)
    print("Columns:")
    print("\n".join(cols))

    wave_col = find_column(cols, ["wavelength", "lambda", "wave"])
    if wave_col is None:
        raise RuntimeError("Could not identify wavelength column. Edit script manually.")
    wave = np.array(tab[wave_col], dtype=float)
    # Guess unit. CIGALE SED wavelengths are often in nm; convert if values look like nm or Angstrom.
    med = np.nanmedian(wave)
    if med > 1000:       # Angstrom likely
        wave_micron = wave / 1e4
    elif med > 100:      # nm likely
        wave_micron = wave / 1000.0
    else:
        wave_micron = wave

    agn_col = find_column(cols, AGN_KEYS)
    total_col = find_column(cols, TOTAL_KEYS, avoid=["stellar", "dust", "agn", "fritz", "skirtor", "torus"])
    if agn_col is None or total_col is None:
        print("Could not confidently identify AGN or total columns.")
        print("Edit AGN_KEYS / TOTAL_KEYS or manually choose columns from the list above.")
        return

    agn = np.array(tab[agn_col], dtype=float)
    total = np.array(tab[total_col], dtype=float)
    print(f"Using wavelength column: {wave_col}")
    print(f"Using AGN column: {agn_col}")
    print(f"Using total column: {total_col}")
    for lam in REST_WAVELENGTHS_MICRON:
        fa = interp_at(wave_micron, agn, lam)
        ft = interp_at(wave_micron, total, lam)
        frac = fa / ft if np.isfinite(fa) and np.isfinite(ft) and ft != 0 else np.nan
        print(f"f_AGN({lam:.2f} micron) = {frac:.3f}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python estimate_optical_agn_fraction.py PATH_TO_best_model.fits")
        raise SystemExit(1)
    main(sys.argv[1])
