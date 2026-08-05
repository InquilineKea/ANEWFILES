#!/usr/bin/env python3
"""
Cell-type frequencies in hippocampus and endothelium, from Allen Brain Cell
Atlas MERFISH (C57BL6J-638850), which has real spatial coordinates.

    python3 abc_fetch.py --dataset merfish          # download (1.6 GB)
    python3 abc_frequencies.py merfish_cells.csv

Produces:
  * frequency tables for hippocampal formation (HPF) by subregion
  * endothelial frequency across brain divisions
  * two .h5ad files carrying REAL obsm['spatial'] coordinates

The statistical caveat that governs everything here
---------------------------------------------------
This dataset is ONE animal (donor C57BL6J-638850), 53 sections. So:

  - Frequencies are descriptive of this brain. They are NOT estimates of a
    population mean, and no confidence interval computed here generalises to
    "mice". The only replication unit available is the section, and sections
    within one brain are pseudoreplicates for any cross-animal claim.
  - Section-level spread is therefore reported as a DISPERSION measure, not as
    a standard error. It tells you how much a frequency wobbles depending on
    where you cut, which is a real and often large effect for a structure with
    strong rostrocaudal gradients.
  - Frequencies are COMPOSITIONAL: they sum to 1 within a region, so they are
    not independent. If one type goes up, others must go down. Do not run
    ordinary correlation or t-tests on raw proportions across regions; use a
    log-ratio transform (CLR/ILR) first. `clr_transform` is provided.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

USECOLS = [
    "cell_label", "brain_section_label", "donor_label", "donor_sex",
    "class", "subclass", "x_reconstructed", "y_reconstructed", "z_reconstructed",
    "parcellation_division", "parcellation_structure",
]
CATCOLS = [
    "class", "subclass", "parcellation_division", "parcellation_structure",
    "brain_section_label", "donor_label", "donor_sex",
]

# Hippocampal formation subregions, in the conventional trisynaptic order.
HPF_ORDER = ["DG", "CA3", "CA2", "CA1", "SUB", "ProS", "POST", "PRE", "PAR",
             "ENTl", "ENTm", "APr", "HATA", "IG", "FC"]

# ABC atlas vascular / perivascular subclasses. "NN" = non-neuronal.
VASCULAR = {
    "333 Endo NN": "endothelial",
    "331 Peri NN": "pericyte",
    "332 SMC NN": "smooth muscle",
    "330 VLMC NN": "vascular leptomeningeal",
    "329 ABC NN": "arachnoid barrier",
}
ENDOTHELIAL = "333 Endo NN"


def load(path: str) -> pd.DataFrame:
    """Load the ABC MERFISH cell metadata, with a pickle cache alongside it."""
    cache = os.path.splitext(path)[0] + "_slim.pkl"
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    df = pd.read_csv(path, usecols=USECOLS, dtype={c: "category" for c in CATCOLS})
    df.to_pickle(cache)
    return df


def clr_transform(props: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """Centred log-ratio transform, for stats on compositional frequencies.

    Proportions live on a simplex, not in Euclidean space. Applying ordinary
    correlations/tests to them produces spurious negative associations purely
    because the parts must sum to 1. CLR maps them somewhere those tools are
    meaningful.
    """
    p = np.asarray(props, dtype=float) + pseudocount
    p = p / p.sum(axis=-1, keepdims=True)
    return np.log(p) - np.log(p).mean(axis=-1, keepdims=True)


def frequency_table(df: pd.DataFrame, group: str, level: str,
                    min_cells: int = 0) -> pd.DataFrame:
    """Counts and within-group proportions of `level` within each `group`."""
    sub = df[[group, level]].astype(str)
    ct = pd.crosstab(sub[group], sub[level])
    if min_cells:
        ct = ct.loc[ct.sum(axis=1) >= min_cells]
    prop = ct.div(ct.sum(axis=1), axis=0)
    return ct, prop


def section_dispersion(df: pd.DataFrame, region_col: str, region: str,
                       level: str, min_cells: int = 200) -> pd.DataFrame:
    """Per-section proportions within one region -> median and IQR.

    This is a DISPERSION measure across sections of a single brain. It is not
    a standard error and does not support cross-animal inference.
    """
    sub = df[df[region_col].astype(str) == region]
    rows = []
    for sec, g in sub.groupby("brain_section_label", observed=True):
        if len(g) < min_cells:
            continue
        vc = g[level].astype(str).value_counts(normalize=True)
        rows.append(vc.rename(sec))
    if not rows:
        return pd.DataFrame()
    mat = pd.DataFrame(rows).fillna(0.0)
    out = pd.DataFrame({
        "median_prop": mat.median(),
        "iqr_lo": mat.quantile(0.25),
        "iqr_hi": mat.quantile(0.75),
        "n_sections": (mat > 0).sum(),
    }).sort_values("median_prop", ascending=False)
    return out


def write_h5ad(df: pd.DataFrame, path: str, note: str) -> None:
    """Write an obs-only AnnData carrying real spatial coordinates.

    There is no expression matrix here -- this is the cell metadata view. X is
    an explicit 0-feature matrix so the object is valid and nobody mistakes a
    placeholder for counts. Fetch the matching expression matrix separately if
    you need genes.
    """
    import anndata as ad
    from scipy import sparse

    obs = df.reset_index(drop=True).copy()
    coords = obs[["x_reconstructed", "y_reconstructed", "z_reconstructed"]].to_numpy(float)
    obs.index = pd.Index(obs["cell_label"].astype(str), name="cell_label")
    adata = ad.AnnData(X=sparse.csr_matrix((len(obs), 0), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = coords[:, :2]
    adata.obsm["spatial_3d"] = coords
    adata.uns["source"] = {
        "atlas": "Allen Brain Cell Atlas",
        "dataset": "MERFISH-C57BL6J-638850-CCF (20231215)",
        "donor": "C57BL6J-638850 (single animal)",
        "coordinates": "x/y/z_reconstructed, millimetres, CCF-registered",
        "note": note,
        "caveat": "Single donor. Section-level spread is dispersion, not a "
                  "standard error; it does not generalise across animals.",
    }
    adata.write_h5ad(path, compression="gzip")
    print("  wrote %s (%d cells)" % (path, adata.n_obs))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("csv", help="cell_metadata_with_parcellation_annotation.csv")
    p.add_argument("--outdir", default=".")
    p.add_argument("--level", default="subclass", choices=["class", "subclass"])
    p.add_argument("--top", type=int, default=12)
    p.add_argument("--no-h5ad", action="store_true")
    args = p.parse_args()

    print("loading %s ..." % args.csv)
    df = load(args.csv)
    print("  %d cells, %d sections, donor(s): %s"
          % (len(df), df.brain_section_label.nunique(),
             ", ".join(map(str, df.donor_label.unique()))))
    if df.donor_label.nunique() == 1:
        print("  NOTE: single donor -> no cross-animal inference is possible here.")

    # ---------------- hippocampus ----------------
    hpf = df[df.parcellation_division.astype(str) == "HPF"].copy()
    print("\n" + "=" * 74)
    print("HIPPOCAMPAL FORMATION (HPF): %d cells" % len(hpf))
    print("=" * 74)

    ct, prop = frequency_table(hpf, "parcellation_structure", args.level)
    order = [r for r in HPF_ORDER if r in prop.index] + \
            [r for r in prop.index if r not in HPF_ORDER]
    prop = prop.loc[order]
    ct = ct.loc[order]

    print("\n--- overall %s frequencies in HPF (top %d) ---" % (args.level, args.top))
    overall = hpf[args.level].astype(str).value_counts(normalize=True)
    overall_n = hpf[args.level].astype(str).value_counts()
    for name, f in overall.head(args.top).items():
        print("  %-34s %7.3f%%  (n=%7d)" % (name, 100 * f, overall_n[name]))

    print("\n--- dominant %s per HPF subregion ---" % args.level)
    print("  %-16s %8s   %-34s %8s" % ("region", "n_cells", "top type", "share"))
    for r in prop.index:
        row = prop.loc[r]
        top = row.idxmax()
        print("  %-16s %8d   %-34s %7.2f%%"
              % (r, int(ct.loc[r].sum()), top, 100 * row[top]))

    # CA1 in detail, with honest section-level dispersion.
    for region in ("CA1", "DG"):
        if region not in prop.index:
            continue
        print("\n--- %s composition, per-section dispersion (top %d) ---"
              % (region, args.top))
        disp = section_dispersion(hpf, "parcellation_structure", region, args.level)
        print("  %-34s %9s %9s %9s %6s"
              % (args.level, "median", "q25", "q75", "nsec"))
        for name, row in disp.head(args.top).iterrows():
            print("  %-34s %8.3f%% %8.3f%% %8.3f%% %6d"
                  % (name, 100 * row.median_prop, 100 * row.iqr_lo,
                     100 * row.iqr_hi, int(row.n_sections)))

    # ---------------- vasculature ----------------
    print("\n" + "=" * 74)
    print("VASCULAR / ENDOTHELIAL")
    print("=" * 74)
    sc = df.subclass.astype(str)
    print("\n--- vascular subclasses, brain-wide ---")
    for key, label in VASCULAR.items():
        n = int((sc == key).sum())
        print("  %-16s %-24s %8d  %6.3f%% of all cells"
              % (key, label, n, 100 * n / len(df)))

    print("\n--- endothelial (%s) share by brain division ---" % ENDOTHELIAL)
    div = df.parcellation_division.astype(str)
    tot = div.value_counts()
    endo = div[sc == ENDOTHELIAL].value_counts()
    share = (endo / tot).dropna().sort_values(ascending=False)
    print("  %-16s %10s %10s %9s" % ("division", "endo", "total", "share"))
    for d, s in share.head(15).items():
        print("  %-16s %10d %10d %8.3f%%" % (d, int(endo.get(d, 0)), int(tot[d]), 100 * s))

    print("\n--- endothelial share within HPF subregions ---")
    hsc = hpf.subclass.astype(str)
    hst = hpf.parcellation_structure.astype(str)
    htot = hst.value_counts()
    hendo = hst[hsc == ENDOTHELIAL].value_counts()
    hshare = (hendo / htot).dropna().sort_values(ascending=False)
    print("  %-16s %9s %9s %9s" % ("region", "endo", "total", "share"))
    for r, s in hshare.items():
        print("  %-16s %9d %9d %8.3f%%" % (r, int(hendo.get(r, 0)), int(htot[r]), 100 * s))

    # ---------------- outputs ----------------
    if not args.no_h5ad:
        print("\nwriting h5ad with real coordinates ...")
        write_h5ad(hpf, os.path.join(args.outdir, "abc_hippocampus.h5ad"),
                   "Hippocampal formation (parcellation_division == HPF)")
        endo_df = df[sc == ENDOTHELIAL].copy()
        write_h5ad(endo_df, os.path.join(args.outdir, "abc_endothelial.h5ad"),
                   "Endothelial cells (subclass == '333 Endo NN'), brain-wide")

    print("\nReminder: frequencies are compositional (they sum to 1 within a "
          "region).\nUse clr_transform() before correlating or testing them.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
