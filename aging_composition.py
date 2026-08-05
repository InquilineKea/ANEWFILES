#!/usr/bin/env python3
"""
Does endothelial / immune cell composition change with age?

Data: Allen Brain Cell Atlas, Zeng-Aging-Mouse-10Xv3. 108 donors in two age
groups (~P56 "adult" vs ~P550 "aged"). Unlike the MERFISH dataset, this has
real biological replication, so cross-animal inference is actually possible --
but it has no spatial coordinates.

    python3 abc_fetch.py --dataset aging-mouse
    python3 aging_composition.py aging_cells.csv

Three corrections without which this analysis returns nonsense
--------------------------------------------------------------
1. FACS. 70 of 108 donors have libraries sorted for `Snap25+ neurons`, which
   depletes exactly the cells being counted here: endothelium 0.51% and
   microglia 0.32% in sorted libraries versus 6.30% and 8.64% in unbiased ones,
   a 12-27x artefact. Only `population_sampling == "unbiased"` is used.

2. Region imbalance. Regions are NOT sampled equally across age groups
   (midbrain is 25.7% of adult cells but 13.1% of aged; isocortex 18.4% vs
   28.2%). Since non-neuronal fractions vary strongly by region, a whole-brain
   adult-vs-aged comparison confounds region with age. Everything here is
   computed WITHIN region and then combined across regions.

3. Compositionality. Proportions within a donor sum to 1, so they are not
   independent -- if neurons are lost, every non-neuronal fraction rises
   mechanically without any cell being gained. Tests run on centred log-ratio
   (CLR) values. The raw percentages are reported too, because they are what
   people want to see, but the inference is on CLR.

The unit of replication is the DONOR, not the cell. Testing 785,769 cells as
though independent would yield absurdly small p-values for trivial effects.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

FOCUS = {
    "333 Endo NN": "Endothelial",
    "331 Peri NN": "Pericyte",
    "332 SMC NN": "Smooth muscle",
    "330 VLMC NN": "VLMC",
    "329 ABC NN": "Arachnoid barrier",
    "334 Microglia NN": "Microglia",
    "335 BAM NN": "Border-assoc. macrophage",
    "338 Lymphoid NN": "Lymphoid",
    "337 DC NN": "Dendritic cell",
}
IMMUNE = ["334 Microglia NN", "335 BAM NN", "338 Lymphoid NN", "337 DC NN"]
VASCULAR = ["333 Endo NN", "331 Peri NN", "332 SMC NN", "330 VLMC NN", "329 ABC NN"]


def clr(props: np.ndarray, pseudocount: float = 1e-4) -> np.ndarray:
    p = np.asarray(props, float) + pseudocount
    p = p / p.sum(axis=-1, keepdims=True)
    lg = np.log(p)
    return lg - lg.mean(axis=-1, keepdims=True)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Non-parametric effect size in [-1, 1]. Robust to the small n here."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if ok.sum() == 0:
        return out
    q = p[ok]
    n = len(q)
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def load(csv: str) -> pd.DataFrame:
    cache = os.path.splitext(csv)[0] + "_analysis.pkl"
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    cols = ["cell_label", "cluster_alias", "donor_label", "donor_age",
            "donor_age_category", "anatomical_division_label", "donor_sex",
            "population_sampling"]
    # cluster_alias must stay integer: it is the merge key against the taxonomy
    # table. Casting it to category makes the merge silently produce all-NaN
    # subclass labels rather than failing, which is why it is spelled out here.
    cat = {c: "category" for c in cols if c not in ("cell_label", "cluster_alias")}
    df = pd.read_csv(csv, usecols=cols, dtype=cat)
    piv = pd.read_csv("cluster_mapping_pivot.csv",
                      usecols=["cluster_alias", "class_name", "subclass_name"])
    df = df.merge(piv, on="cluster_alias", how="left")
    # Fail loudly on a broken join. A dtype mismatch on the merge key yields
    # all-NaN labels and an empty, silently wrong analysis downstream.
    n_bad = int(df.subclass_name.isna().sum())
    if n_bad:
        raise RuntimeError(
            "taxonomy join failed for %d/%d cells (subclass_name is NaN). "
            "Check that cluster_alias dtypes match on both sides."
            % (n_bad, len(df)))
    df = df[df.population_sampling == "unbiased"].copy()
    df.to_pickle(cache)
    return df


def donor_region_props(df: pd.DataFrame, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Proportion of each subclass within each donor x region."""
    g = df.groupby(["donor_label", "anatomical_division_label", "donor_age_category"],
                   observed=True)
    rows, meta = [], []
    for (donor, region, age), sub in g:
        if len(sub) < min_cells:
            continue
        vc = sub.subclass_name.astype(str).value_counts(normalize=True)
        rows.append(vc)
        meta.append({"donor": donor, "region": region, "age": age, "n_cells": len(sub)})
    props = pd.DataFrame(rows).fillna(0.0).reset_index(drop=True)
    return props, pd.DataFrame(meta)


def test_within_region(props: pd.DataFrame, meta: pd.DataFrame,
                       targets: list[str], min_donors: int = 4) -> pd.DataFrame:
    """Adult vs aged, per subclass per region, on CLR values."""
    clr_all = pd.DataFrame(clr(props.to_numpy()), columns=props.columns)
    out = []
    for region in meta.region.unique():
        m = (meta.region == region).to_numpy()
        ad = m & (meta.age == "adult").to_numpy()
        ag = m & (meta.age == "aged").to_numpy()
        if ad.sum() < min_donors or ag.sum() < min_donors:
            continue
        for t in targets:
            if t not in props.columns:
                continue
            a, b = clr_all.loc[ad, t].to_numpy(), clr_all.loc[ag, t].to_numpy()
            pa, pb = props.loc[ad, t].to_numpy(), props.loc[ag, t].to_numpy()
            try:
                u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
            except ValueError:
                continue
            out.append({
                "region": region, "subclass": t, "label": FOCUS.get(t, t),
                "n_adult": int(ad.sum()), "n_aged": int(ag.sum()),
                "pct_adult": 100 * np.median(pa), "pct_aged": 100 * np.median(pb),
                "log2fc": np.log2((np.median(pb) + 1e-6) / (np.median(pa) + 1e-6)),
                "clr_diff": float(np.median(b) - np.median(a)),
                "cliffs_delta": cliffs_delta(b, a), "p": p,
            })
    res = pd.DataFrame(out)
    if len(res):
        res["q"] = bh_fdr(res.p.to_numpy())
    return res


def test_pooled(props: pd.DataFrame, meta: pd.DataFrame,
                targets: list[str]) -> pd.DataFrame:
    """Region-stratified overall test.

    Each donor x region observation is centred on its region's adult median
    before pooling, so region differences cannot masquerade as age effects.
    """
    clr_all = pd.DataFrame(clr(props.to_numpy()), columns=props.columns)
    out = []
    for t in targets:
        if t not in props.columns:
            continue
        adj_a, adj_b = [], []
        for region in meta.region.unique():
            m = (meta.region == region).to_numpy()
            ad = m & (meta.age == "adult").to_numpy()
            ag = m & (meta.age == "aged").to_numpy()
            if ad.sum() < 3 or ag.sum() < 3:
                continue
            base = np.median(clr_all.loc[ad, t])
            adj_a.append(clr_all.loc[ad, t].to_numpy() - base)
            adj_b.append(clr_all.loc[ag, t].to_numpy() - base)
        if not adj_a:
            continue
        a, b = np.concatenate(adj_a), np.concatenate(adj_b)
        u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
        out.append({
            "subclass": t, "label": FOCUS.get(t, t),
            "n_adult_obs": len(a), "n_aged_obs": len(b),
            "clr_shift": float(np.median(b) - np.median(a)),
            "cliffs_delta": cliffs_delta(b, a), "p": p,
        })
    res = pd.DataFrame(out)
    if len(res):
        res["q"] = bh_fdr(res.p.to_numpy())
    return res.sort_values("clr_shift", ascending=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="aging_cells.csv")
    ap.add_argument("--min-cells", type=int, default=150,
                    help="min cells per donor x region to include that observation")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    df = load(args.csv)
    print("unbiased cells: %d | donors: %d (adult %d, aged %d)"
          % (len(df), df.donor_label.nunique(),
             df[df.donor_age_category == "adult"].donor_label.nunique(),
             df[df.donor_age_category == "aged"].donor_label.nunique()))

    props, meta = donor_region_props(df, args.min_cells)
    print("donor x region observations: %d (>= %d cells each)" % (len(meta), args.min_cells))

    targets = [t for t in FOCUS if t in props.columns]

    print("\n" + "=" * 78)
    print("REGION-STRATIFIED OVERALL TEST (unit = donor x region, CLR)")
    print("=" * 78)
    pooled = test_pooled(props, meta, targets)
    print("  %-26s %8s %8s %9s %10s %6s"
          % ("cell type", "CLRshift", "delta", "p", "q(FDR)", "sig"))
    for _, r in pooled.iterrows():
        sig = "***" if r.q < 0.001 else "**" if r.q < 0.01 else "*" if r.q < 0.05 else "ns"
        print("  %-26s %+8.3f %+8.3f %9.2e %10.2e %6s"
              % (r.label, r.clr_shift, r.cliffs_delta, r.p, r.q, sig))

    print("\n" + "=" * 78)
    print("PER-REGION TESTS (aged vs adult)")
    print("=" * 78)
    per = test_within_region(props, meta, targets)
    sig = per[per.q < 0.05].sort_values("q")
    print("  %d of %d region x celltype tests significant at FDR 5%%" % (len(sig), len(per)))
    if len(sig):
        print("\n  %-10s %-26s %8s %8s %8s %9s"
              % ("region", "cell type", "%adult", "%aged", "delta", "q"))
        for _, r in sig.iterrows():
            print("  %-10s %-26s %7.2f%% %7.2f%% %+8.3f %9.2e"
                  % (r.region, r.label, r.pct_adult, r.pct_aged, r.cliffs_delta, r.q))

    print("\n" + "=" * 78)
    print("IMMUNE COMPOSITION WITHIN THE IMMUNE COMPARTMENT")
    print("=" * 78)
    imm = df[df.class_name.astype(str).str.contains("Immune", na=False)]
    ip = (imm.groupby(["donor_label", "donor_age_category"], observed=True)
            .subclass_name.value_counts(normalize=True).unstack(fill_value=0.0))
    ip = ip[[c for c in IMMUNE if c in ip.columns]]
    idx = ip.index.get_level_values("donor_age_category")
    print("  %-26s %10s %10s %9s %9s"
          % ("immune subtype", "%adult", "%aged", "p", "q"))
    ps, rows = [], []
    for c in ip.columns:
        a, b = ip.loc[idx == "adult", c].to_numpy(), ip.loc[idx == "aged", c].to_numpy()
        u, p = stats.mannwhitneyu(b, a, alternative="two-sided")
        ps.append(p)
        rows.append((FOCUS.get(c, c), 100 * np.median(a), 100 * np.median(b), p))
    qs = bh_fdr(np.array(ps))
    for (lab, ma, mb, p), q in zip(rows, qs):
        star = "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"
        print("  %-26s %9.3f%% %9.3f%% %9.2e %9.2e %s" % (lab, ma, mb, p, q, star))

    props.assign(**meta).to_csv(os.path.join(args.outdir, "aging_donor_region_props.csv"),
                                index=False)
    pooled.to_csv(os.path.join(args.outdir, "results_aging_pooled.csv"), index=False)
    per.to_csv(os.path.join(args.outdir, "results_aging_per_region.csv"), index=False)
    ip.reset_index().to_csv(os.path.join(args.outdir, "aging_immune_props.csv"), index=False)
    print("\nwrote result CSVs to %s" % args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
