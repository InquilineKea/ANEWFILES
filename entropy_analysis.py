#!/usr/bin/env python3
"""
Does spatial entropy increase with age?

    python3 entropy_analysis.py

Short answer: **that exact question cannot be tested with any dataset reachable
here**, and this script is explicit about why rather than substituting a
lookalike and calling it an answer.

Spatial entropy -- how mixed cell types are in a local neighbourhood -- needs
per-cell coordinates. Every spatially-resolved dataset in the Allen Brain Cell
Atlas is a SINGLE donor, and they are different species:

    MERFISH C57BL6J-638850   mouse     1 donor
    Xenium CJ23.56.004       marmoset  1 donor, 5.7 yr
    MERSCOPE H22.30.001      human     1 donor, 50 yr
    MERSCOPE QM23.50.001     macaque   1 donor, 10 yr

Contrasting those would confound species with age completely. Meanwhile the
dataset with 108 donors and a real age contrast (Zeng-Aging-Mouse-10Xv3) is
dissociated scRNA-seq: the coordinates do not exist.

So this script does the three things that ARE legitimate:

  A. Builds and validates a neighbourhood spatial-entropy measure on the real
     MERFISH coordinates, checking it recovers known anatomy (it must score the
     dentate granule layer as highly ordered and neuropil-rich regions as
     mixed). This is the instrument, ready for data that has both axes.

  B. Tests COMPOSITIONAL entropy against age in the 108-donor dataset. This is
     a genuinely different quantity -- diversity of cell types in a sample,
     with no geometry -- and is reported as such. It is the closest testable
     relative of the question, not a substitute for it.

  C. Reports what a real test would need, via the effect size B could detect.

Two estimator traps handled
---------------------------
* Shannon entropy is strongly biased by sample size: more cells sampled ->
  higher observed entropy. Aged and adult samples do not contain equal numbers
  of cells, so raw entropy would partly measure sequencing depth. Every sample
  is RAREFIED to an identical cell count, averaged over repeated draws.
* The neuronal-dropout confound from `aging_confound_check.py` applies here
  too: if aged samples are relatively depleted of neurons, the cell-type mix
  shifts and entropy moves with it for non-biological reasons. Entropy is
  therefore computed twice -- over all cells, and within the non-neuronal
  compartment alone.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aging_composition import bh_fdr, cliffs_delta, load  # noqa: E402
from aging_confound_check import is_non_neuronal  # noqa: E402

MERFISH_CSV = "merfish_cells.csv"
AGING_CSV = "aging_cells.csv"


def shannon(counts: np.ndarray) -> float:
    """Shannon entropy in nats of a count vector."""
    c = np.asarray(counts, float)
    c = c[c > 0]
    if c.sum() <= 0:
        return np.nan
    p = c / c.sum()
    return float(-(p * np.log(p)).sum())


def rarefied_entropy(labels: np.ndarray, n: int, reps: int = 30,
                     seed: int = 0) -> float:
    """Mean Shannon entropy over `reps` draws of exactly `n` labels.

    Rarefaction is not optional here. Observed entropy rises with the number of
    cells sampled, so comparing samples of different size measures depth as
    much as diversity.
    """
    if len(labels) < n:
        return np.nan
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(reps):
        pick = rng.choice(labels, size=n, replace=False)
        _, counts = np.unique(pick, return_counts=True)
        vals.append(shannon(counts))
    return float(np.mean(vals))


# ------------------------------------------------------------------ part A
def spatial_entropy_merfish(k: int = 20, max_cells_per_region: int = 20000) -> pd.DataFrame:
    """Neighbourhood cell-type entropy on real MERFISH coordinates.

    For each cell: take its k nearest neighbours WITHIN THE SAME SECTION
    (sections are separate physical slices; mixing them would create neighbours
    that are millimetres apart in z), then take the Shannon entropy of the
    subclass composition of that neighbourhood.

    Low entropy = locally homogeneous (a tight, pure layer).
    High entropy = locally mixed.
    """
    print("[A] spatial entropy on MERFISH coordinates (k=%d neighbours)" % k)
    cols = ["brain_section_label", "subclass", "parcellation_division",
            "parcellation_structure", "x_reconstructed", "y_reconstructed"]
    df = pd.read_csv(MERFISH_CSV, usecols=cols,
                     dtype={c: "category" for c in cols if not c.startswith(("x_", "y_"))})
    print("    %d cells, %d sections" % (len(df), df.brain_section_label.nunique()))

    rng = np.random.default_rng(0)
    out = []
    for sec, g in df.groupby("brain_section_label", observed=True):
        if len(g) < k + 1:
            continue
        xy = g[["x_reconstructed", "y_reconstructed"]].to_numpy(float)
        codes = g.subclass.astype(str).to_numpy()
        tree = cKDTree(xy)
        # neighbour indices, excluding self (column 0)
        _, idx = tree.query(xy, k=k + 1)
        neigh = codes[idx[:, 1:]]
        ent = np.array([shannon(np.unique(row, return_counts=True)[1]) for row in neigh])
        out.append(pd.DataFrame({
            "section": sec, "division": g.parcellation_division.astype(str).to_numpy(),
            "structure": g.parcellation_structure.astype(str).to_numpy(),
            "subclass": codes, "entropy": ent}))
    E = pd.concat(out, ignore_index=True)
    print("    computed neighbourhood entropy for %d cells" % len(E))
    return E


# ------------------------------------------------------------------ part B
def compositional_entropy_aging(n_rarefy: int = 300, reps: int = 30) -> pd.DataFrame:
    """Cell-type diversity per donor x region, adult vs aged. NOT spatial."""
    print("\n[B] compositional entropy vs age (108-donor dataset, no geometry)")
    df = load(AGING_CSV)
    nn = is_non_neuronal(df.class_name)
    rows = []
    for scope, sub in (("all cells", df), ("non-neuronal only", df[nn])):
        for (donor, region, age), s in sub.groupby(
                ["donor_label", "anatomical_division_label", "donor_age_category"],
                observed=True):
            lab = s.subclass_name.astype(str).to_numpy()
            e = rarefied_entropy(lab, n_rarefy, reps)
            if np.isnan(e):
                continue
            rows.append({"scope": scope, "donor": donor, "region": region,
                         "age": age, "entropy": e, "n_cells": len(s),
                         "n_types": len(np.unique(lab))})
    return pd.DataFrame(rows)


def test_entropy_vs_age(E: pd.DataFrame) -> pd.DataFrame:
    """Region-stratified adult vs aged test on entropy."""
    out = []
    for scope, sub in E.groupby("scope"):
        A, B = [], []
        for region, g in sub.groupby("region", observed=True):
            a = g[g.age == "adult"].entropy.to_numpy()
            b = g[g.age == "aged"].entropy.to_numpy()
            if len(a) < 3 or len(b) < 3:
                continue
            base = np.median(a)
            A.append(a - base)
            B.append(b - base)
        if not A:
            continue
        A, B = np.concatenate(A), np.concatenate(B)
        p = stats.mannwhitneyu(B, A, alternative="two-sided")[1]
        out.append({"scope": scope, "n_adult": len(A), "n_aged": len(B),
                    "median_adult": float(np.median(sub[sub.age == "adult"].entropy)),
                    "median_aged": float(np.median(sub[sub.age == "aged"].entropy)),
                    "shift_nats": float(np.median(B) - np.median(A)),
                    "delta": cliffs_delta(B, A), "p": p})
    R = pd.DataFrame(out)
    if len(R):
        R["q"] = bh_fdr(R.p.to_numpy())
    return R


def power_note(E: pd.DataFrame, scope: str = "non-neuronal only") -> None:
    """What effect size would this design actually detect?"""
    sub = E[E.scope == scope]
    sd = float(sub.entropy.std())
    n_a = int((sub.age == "adult").sum())
    n_b = int((sub.age == "aged").sum())
    # Two-sample t-test detectable difference at 80% power, alpha 0.05.
    d80 = 2.8 * sd * np.sqrt(1 / max(n_a, 1) + 1 / max(n_b, 1))
    print("\n[C] sensitivity (%s)" % scope)
    print("    entropy SD across samples: %.3f nats" % sd)
    print("    n = %d adult / %d aged observations" % (n_a, n_b))
    print("    smallest difference detectable at ~80%% power: %.3f nats "
          "(%.1f%% of SD)" % (d80, 100 * d80 / sd))


def main() -> int:
    if not os.path.exists(AGING_CSV):
        print("missing %s -- run: python3 abc_fetch.py --dataset aging-mouse" % AGING_CSV)
        return 2

    # ---- B: the testable age question ----
    E = compositional_entropy_aging()
    R = test_entropy_vs_age(E)
    print("\n    compositional entropy, aged vs adult (region-stratified, rarefied):")
    print("    %-20s %9s %9s %9s %8s %9s %s"
          % ("scope", "adult", "aged", "shift", "delta", "q", ""))
    for _, r in R.iterrows():
        sig = "***" if r.q < 0.001 else "**" if r.q < 0.01 else "*" if r.q < 0.05 else "ns"
        print("    %-20s %9.3f %9.3f %+9.3f %+8.3f %9.2e %s"
              % (r.scope, r.median_adult, r.median_aged, r.shift_nats,
                 r.delta, r.q, sig))
    E.to_csv("results_entropy_aging.csv", index=False)
    R.to_csv("results_entropy_aging_test.csv", index=False)
    power_note(E)

    # ---- A: the instrument, on real coordinates ----
    if os.path.exists(MERFISH_CSV):
        S = spatial_entropy_merfish()
        by_div = (S.groupby("division").entropy
                  .agg(["mean", "median", "count"])
                  .sort_values("median", ascending=False))
        print("\n    neighbourhood entropy by brain division (nats, k=20):")
        for d, r in by_div.head(8).iterrows():
            print("      %-16s %.3f  (n=%d)" % (d, r["median"], r["count"]))
        print("      ...")
        for d, r in by_div.tail(3).iterrows():
            print("      %-16s %.3f  (n=%d)" % (d, r["median"], r["count"]))

        hpf = S[S.division == "HPF"]
        print("\n    within hippocampal formation, by subregion:")
        for st, r in (hpf.groupby("structure").entropy
                      .agg(["median", "count"]).sort_values("median").iterrows()):
            if r["count"] < 500:
                continue
            print("      %-16s %.3f  (n=%d)" % (st, r["median"], r["count"]))
        S.to_csv("results_entropy_spatial_merfish.csv.gz", index=False,
                 compression="gzip")
        print("\n    wrote results_entropy_spatial_merfish.csv.gz")
    else:
        print("\n[A] skipped: %s not present" % MERFISH_CSV)
    return 0


if __name__ == "__main__":
    sys.exit(main())
