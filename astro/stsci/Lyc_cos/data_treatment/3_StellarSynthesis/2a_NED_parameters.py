#!/usr/bin/env python
"""
ned_distances_extinction.py

Retrieve distances (Hubble flow and redshift-independent) and foreground
Galactic extinction from NED for a list of galaxies.

Uses NED's current REST API (https://ned.ipac.caltech.edu/Docs::API/).
The legacy cgi-bin endpoints that astroquery.ipac.ned targets were deprecated
in January 2026, so this talks to the new services directly. Every response is
a VOTable, so no HTML scraping is needed.

Just edit GALAXY_NAMES below and run the script (or import build_catalogue
and call it yourself) -- no command-line arguments needed.
"""

import time
import warnings
from io import BytesIO

import numpy as np
import requests
from astropy.io.votable import parse as parse_votable
from astropy.io.votable.exceptions import VOWarning
from astropy.table import Table, vstack
from pathlib import Path
import lime

NED_API = "https://ned.ipac.caltech.edu/NED::API"
SESSION = requests.Session()


# ----------------------------------------------------------------------------
# EDIT THIS LIST
# ----------------------------------------------------------------------------
GALAXY_NAMES = [
    'SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
    'UGCA281', 'NGC1705', 'Haro2', 'NGC4861', 'He2-10', 'MRK1450',
    'UGC4483', 'VIIZw403', 'NGC2366',
    'Haro11_A', 'Haro11_B', 'Haro11_C',
    'IZw18_SE']

# NED resolves catalog/IAU names, not sub-region labels. Entries here are
# rewritten to the NED-resolvable name before querying; the _A/_B/_C suffix
# is kept in the output's input_name column so you can still tell them apart.
NAME_OVERRIDES = {'Haro11_A': 'Haro 11',
                  'Haro11_B': 'Haro 11',
                  'Haro11_C': 'Haro 11',
                  'SBS1159545': 'SBS 1159+545',
                  'SBS1415437': 'SBS 1415+437',
                  'SBS0335052': 'SBS 0335-052',
                  'IZw18_SE': 'IZw18'}


OUT_PREFIX = "ned_output"      # -> ned_output.csv, ned_output_distances.csv, ...
WITH_DISTANCES = True          # dump all individual redshift-indep. estimates
WITH_BANDS = True              # dump per-bandpass extinction
PAUSE = 1.0                    # seconds between queries (be polite to NED)


# ----------------------------------------------------------------------------
# low-level access
# ----------------------------------------------------------------------------
def ned_query(service, timeout=60, n_retries=3, backoff=2.0, **params):
    """Call a NED API service; return every TABLE in the VOTable as a list."""
    url = f"{NED_API}/{service}"
    for attempt in range(n_retries):
        try:
            resp = SESSION.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            break
        except requests.RequestException as exc:
            if attempt == n_retries - 1:
                raise RuntimeError(f"{service} failed for {params}: {exc}") from exc
            time.sleep(backoff ** attempt)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", VOWarning)
        vot = parse_votable(BytesIO(resp.content), verify="warn")
        return [t.to_table(use_names_over_ids=True) for t in vot.iter_tables()]


def _scalar(value):
    """Masked / empty -> nan; otherwise a plain python or numpy scalar."""
    if value is None or value is np.ma.masked:
        return np.nan
    if isinstance(value, bytes):
        value = value.decode()
    if isinstance(value, str) and not value.strip():
        return np.nan
    return value


def _flatten(tables):
    """Merge the first row of every returned table into one dict."""
    out = {}
    for tab in tables:
        if len(tab) == 0:
            continue
        for col in tab.colnames:
            out.setdefault(col, _scalar(tab[col][0]))
    return out


def _get(row, *keys, default=np.nan):
    """Case-insensitive exact-then-substring lookup, so minor label changes
    in the NED schema don't break the script."""
    low = {k.lower(): v for k, v in row.items()}
    for key in keys:
        if key.lower() in low:
            return low[key.lower()]
    for key in keys:
        for name, val in low.items():
            if key.lower() in name:
                return val
    return default


# ----------------------------------------------------------------------------
# per-object products
# ----------------------------------------------------------------------------
def object_overview(query_name):
    """One summary row for a single NED-resolvable name."""
    row = _flatten(ned_query("OverviewOfObject", TARGET=query_name))
    return {
        "ned_name": _get(row, "Object Name", "Crossidentifications", default=""),
        "ra": _get(row, "Lon (Equatorial J2000)", "Lon (Equatorial"),
        "dec": _get(row, "Lat (Equatorial J2000)", "Lat (Equatorial"),
        "z": _get(row, "Redshift"),
        "z_err": _get(row, "Redshift Unc"),
        "v_helio": _get(row, "v (Heliocentric)"),
        "D_3K": _get(row, "D (3K CMB)"),          # Mpc, Hubble flow
        "D_3K_err": _get(row, "D Unc (3K CMB)"),
        "D_mean": _get(row, "Mean Distance"),      # Mpc, redshift-independent
        "D_mean_sem": _get(row, "SEM Distance"),
        "A_V": _get(row, "Galactic Extinction (Landolt V)", "Landolt V"),
        "A_K": _get(row, "Galactic Extinction (UKIRT K)", "UKIRT K"),
    }


def redshift_independent_distances(query_name):
    """Full NED-D style table: (m-M), err, D[Mpc], method, refcode, ..."""
    tables = [t for t in ned_query("DistancesOfObject", TARGET=query_name) if len(t)]
    if not tables:
        return None
    return vstack(tables, metadata_conflicts="silent") if len(tables) > 1 else tables[0]


def extinction_bandpasses(query_name):
    """Per-filter foreground Galactic extinction, with central wavelength in um."""
    tables = [t for t in ned_query("ExtinctionAtTarget", TARGET=query_name) if len(t)]
    if not tables:
        return None
    return vstack(tables, metadata_conflicts="silent") if len(tables) > 1 else tables[0]


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def build_catalogue(names, overrides=None, pause=1.0, with_distances=False,
                    with_bands=False, verbose=True):
    """
    names:      list of galaxy identifiers, exactly as you want them labeled
                in the output (e.g. 'Haro11_A').
    overrides:  optional dict mapping an entry in `names` to the actual
                string sent to NED (e.g. {'Haro11_A': 'Haro 11'}); useful
                for sub-region names NED doesn't resolve on its own.
    Returns (summary, distances, bands) astropy Tables. distances/bands are
    None if with_distances/with_bands is False or nothing came back.
    """
    overrides = overrides or {}
    rows, dist_tabs, band_tabs = [], [], []

    for i, name in enumerate(names):
        query_name = overrides.get(name, name)
        try:
            overview = object_overview(query_name)
            rows.append({"input_name": name, "query_name": query_name, **overview})

            if with_distances:
                t = redshift_independent_distances(query_name)
                if t is not None:
                    t["input_name"] = name
                    dist_tabs.append(t)

            if with_bands:
                t = extinction_bandpasses(query_name)
                if t is not None:
                    t["input_name"] = name
                    band_tabs.append(t)

            status = "ok"
        except Exception as exc:
            blank = {k: np.nan for k in (
                "ra", "dec", "z", "z_err", "v_helio", "D_3K", "D_3K_err",
                "D_mean", "D_mean_sem", "A_V", "A_K")}
            rows.append({"input_name": name, "query_name": query_name,
                        "ned_name": "", **blank})
            status = f"FAILED ({exc})"

        if verbose:
            print(f"[{i + 1}/{len(names)}] {name} -> {query_name}: {status}", flush=True)
        if pause:
            time.sleep(pause)

    summary = Table(rows)
    distances = vstack(dist_tabs, metadata_conflicts="silent") if dist_tabs else None
    bands = vstack(band_tabs, metadata_conflicts="silent") if band_tabs else None
    return summary, distances, bands


if __name__ == "__main__":
    summary, distances, bands = build_catalogue(
        GALAXY_NAMES,
        overrides=NAME_OVERRIDES,
        pause=PAUSE,
        with_distances=WITH_DISTANCES,
        with_bands=WITH_BANDS,
    )

    summary.write(f"{OUT_PREFIX}.csv", format="csv", overwrite=True)
    print(f"\nwrote {OUT_PREFIX}.csv")

    if distances is not None:
        distances.write(f"{OUT_PREFIX}_distances.csv", format="csv", overwrite=True)
        print(f"wrote {OUT_PREFIX}_distances.csv")

    if bands is not None:
        bands.write(f"{OUT_PREFIX}_extinction.csv", format="csv", overwrite=True)
        print(f"wrote {OUT_PREFIX}_extinction.csv")