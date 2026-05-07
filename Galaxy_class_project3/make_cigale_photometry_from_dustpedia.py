"""
Build a small CIGALE input table from DustPedia aperture-matched photometry.

This script is intentionally verbose because DustPedia/VizieR column names can
change slightly between access methods. It queries the DustPedia VizieR table
J/A+A/609/A37/apphot and writes a CIGALE-compatible CSV with fluxes in mJy.

Requirements:
    pip install pyvo pandas astropy numpy

Output:
    cigale_dustpedia_input.csv

If the automatic column matching fails, run the script once, inspect the printed
column list, and update BAND_CANDIDATES below.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
import pyvo
from astropy.table import Table

TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
APPHOT_TABLE = '"J/A+A/609/A37/apphot"'
SAMPLE_TABLE = '"J/A+A/609/A37/sample"'

TARGETS = ["NGC0628", "NGC1097", "NGC1068", "NGC4151", "NGC5055"]

# CIGALE filter names used in Boquien et al. 2019 examples and common CIGALE installs.
# If your CIGALE version uses different names, change these labels after running
# `pcigale-filters list`.
CIGALE_BANDS = {
    "galex.FUV": ["FUV", "GALEX_FUV", "galex_fuv", "FUVflux", "S_FUV", "f_FUV"],
    "galex.NUV": ["NUV", "GALEX_NUV", "galex_nuv", "NUVflux", "S_NUV", "f_NUV"],
    "sdss.up": ["u", "umag", "SDSS_u", "sdss_u", "u_SDSS", "F_u", "S_u"],
    "sdss.gp": ["g", "gmag", "SDSS_g", "sdss_g", "g_SDSS", "F_g", "S_g"],
    "sdss.rp": ["r", "rmag", "SDSS_r", "sdss_r", "r_SDSS", "F_r", "S_r"],
    "sdss.ip": ["i", "imag", "SDSS_i", "sdss_i", "i_SDSS", "F_i", "S_i"],
    "sdss.zp": ["z", "zmag", "SDSS_z", "sdss_z", "z_SDSS", "F_z", "S_z"],
    "2mass.J": ["J", "2MASS_J", "2mass_J", "J_2MASS", "F_J", "S_J"],
    "2mass.H": ["H", "2MASS_H", "2mass_H", "H_2MASS", "F_H", "S_H"],
    "2mass.Ks": ["Ks", "K", "2MASS_Ks", "2mass_Ks", "Ks_2MASS", "F_Ks", "S_Ks"],
    "WISE1": ["W1", "WISE1", "WISE_W1", "wise_w1", "F_W1", "S_W1"],
    "WISE2": ["W2", "WISE2", "WISE_W2", "wise_w2", "F_W2", "S_W2"],
    "WISE3": ["W3", "WISE3", "WISE_W3", "wise_w3", "F_W3", "S_W3"],
    "WISE4": ["W4", "WISE4", "WISE_W4", "wise_w4", "F_W4", "S_W4"],
}

# Try to find matching uncertainty columns. If none are found, use a conservative
# fractional uncertainty floor.
ERR_CANDIDATE_SUFFIXES = ["_err", "err", "_e", "e_", "_unc", "unc", "_Error", "_error"]
DEFAULT_FRAC_ERROR = 0.10
MIN_RELATIVE_ERROR = 0.05

@dataclass
class TargetMeta:
    id: str
    redshift: float
    distance_mpc: float

TARGET_META = {
    "NGC0628": TargetMeta("NGC0628", 0.00219, 9.84),
    "NGC1097": TargetMeta("NGC1097", 0.00424, 14.5),
    "NGC1068": TargetMeta("NGC1068", 0.00379, 14.4),
    "NGC4151": TargetMeta("NGC4151", 0.00332, 15.8),
    "NGC5055": TargetMeta("NGC0628", 0.001668, 8.870),
}


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def find_col(columns: Iterable[str], candidates: list[str]) -> str | None:
    by_norm = {norm(c): c for c in columns}
    for cand in candidates:
        if norm(cand) in by_norm:
            return by_norm[norm(cand)]
    # relaxed contains match
    for cand in candidates:
        nc = norm(cand)
        for col in columns:
            if nc and nc in norm(col):
                return col
    return None


def find_err_col(columns: Iterable[str], flux_col: str) -> str | None:
    nflux = norm(flux_col)
    by_norm = {norm(c): c for c in columns}
    guesses = []
    for suffix in ERR_CANDIDATE_SUFFIXES:
        guesses.append(flux_col + suffix)
        guesses.append("e_" + flux_col)
        guesses.append("err_" + flux_col)
    for guess in guesses:
        if norm(guess) in by_norm:
            return by_norm[norm(guess)]
    for col in columns:
        nc = norm(col)
        if nflux in nc and any(token in nc for token in ["err", "unc", "error"]):
            return col
    return None


def query_table() -> pd.DataFrame:
    svc = pyvo.dal.TAPService(TAP_URL)
    # Get a broad row set. We keep the query simple because table metadata can vary.
    query = f"SELECT * FROM {APPHOT_TABLE}"
    print("Querying DustPedia aperture photometry from VizieR TAP...")
    result = svc.search(query)
    tab: Table = result.to_table()
    df = tab.to_pandas()
    print(f"Downloaded {len(df)} rows and {len(df.columns)} columns.")
    print("Columns available in apphot table:")
    print(", ".join(map(str, df.columns)))
    return df


def identify_name_column(df: pd.DataFrame) -> str:
    candidates = ["Name", "name", "Galaxy", "galaxy", "objname", "Object", "PGC", "ID"]
    col = find_col(df.columns, candidates)
    if col is None:
        raise RuntimeError(
            "Could not identify the galaxy-name column. Inspect the printed columns and edit identify_name_column()."
        )
    return col


def select_targets(df: pd.DataFrame, name_col: str) -> pd.DataFrame:
    def compact_name(x):
        return norm(x).replace("ngc", "ngc")

    wanted = {norm(t): t for t in TARGETS}
    rows = []
    for _, row in df.iterrows():
        n = norm(row[name_col])
        if n in wanted:
            rows.append(row)
    out = pd.DataFrame(rows)
    if len(out) != len(TARGETS):
        print("WARNING: did not find all targets by exact normalized name.")
        print("Found:", list(out[name_col]) if len(out) else [])
        print("Trying relaxed contains matching...")
        rows = []
        for target in TARGETS:
            tnorm = norm(target)
            matches = df[df[name_col].astype(str).map(lambda x: tnorm in norm(x))]
            if len(matches) > 0:
                rows.append(matches.iloc[0])
        out = pd.DataFrame(rows)
    if len(out) == 0:
        raise RuntimeError("No target galaxies matched. Inspect name column and TARGETS.")
    return out


def build_cigale_table(df: pd.DataFrame, name_col: str) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        raw_name = str(r[name_col])
        # Normalize to our target identifiers.
        target_id = None
        for t in TARGETS:
            if norm(t) in norm(raw_name):
                target_id = t
                break
        if target_id is None:
            target_id = raw_name.replace(" ", "")
        meta = TARGET_META.get(target_id, TargetMeta(target_id, np.nan, np.nan))
        row = {"id": target_id, "redshift": meta.redshift, "distance": meta.distance_mpc}
        for cigale_name, candidates in CIGALE_BANDS.items():
            flux_col = find_col(df.columns, candidates)
            if flux_col is None:
                row[cigale_name] = np.nan
                row[cigale_name + "_err"] = np.nan
                continue
            flux = r[flux_col]
            try:
                flux = float(flux)
            except Exception:
                flux = np.nan
            # DustPedia VizieR flux columns are expected to be Jy in many tables.
            # If a column unit is already mJy in your VizieR export, change UNIT_SCALE_TO_MJY below.
            # The script writes a metadata note and asks the user to verify units.
            UNIT_SCALE_TO_MJY = 1000.0
            flux_mjy = flux * UNIT_SCALE_TO_MJY if np.isfinite(flux) else np.nan
            err_col = find_err_col(df.columns, flux_col)
            if err_col is not None:
                try:
                    err_mjy = float(r[err_col]) * UNIT_SCALE_TO_MJY
                except Exception:
                    err_mjy = np.nan
            else:
                err_mjy = np.nan
            if not np.isfinite(err_mjy) and np.isfinite(flux_mjy):
                err_mjy = DEFAULT_FRAC_ERROR * abs(flux_mjy)
            if np.isfinite(flux_mjy) and np.isfinite(err_mjy):
                err_mjy = max(err_mjy, MIN_RELATIVE_ERROR * abs(flux_mjy))
            row[cigale_name] = flux_mjy
            row[cigale_name + "_err"] = err_mjy
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    df = query_table()
    name_col = identify_name_column(df)
    print(f"Using name column: {name_col}")
    targets = select_targets(df, name_col)
    cig = build_cigale_table(targets, name_col)
    # Put columns in a clean order.
    ordered = ["id", "redshift", "distance"]
    for b in CIGALE_BANDS:
        ordered += [b, b + "_err"]
    cig = cig[ordered]
    out = "cigale_dustpedia_input.csv"
    cig.to_csv(out, index=False, na_rep="")
    print(f"Wrote {out}")
    print(cig)
    print("\nIMPORTANT: Verify the DustPedia flux units in the VizieR table metadata.")
    print("This script assumes flux columns are Jy and converts to mJy for CIGALE.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERROR:", exc, file=sys.stderr)
        print("Inspect the printed column names and edit the candidate lists near the top of the script.", file=sys.stderr)
        raise
