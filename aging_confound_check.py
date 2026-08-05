#!/usr/bin/env python3
"""
The neuronal-denominator confound, and what survives it.

    python3 aging_confound_check.py aging_cells.csv

Why this file exists
--------------------
`aging_composition.py` computes each cell type's share of ALL cells in a
donor x region sample. That is the obvious thing to do and it is confounded.

Neurons are fragile in dissociation, and their recovery differs by age group:
32.2% of unbiased-library cells are neuronal in adults versus 24.7% in aged
mice (a 23% relative drop). Median neuronal fraction per donor x region falls
much further, 34.1% to 9.5%, because aged samples are more often
non-neuronal-dominated. Library quality does not explain it -- median genes per
neuron is 7548 adult vs 7405 aged, essentially identical -- so this is a change
in which cells survive to be sequenced, not in how well they were measured.

Whether that neuron loss is technical or partly biological cannot be settled
from this data. It does not matter for the correction: proportions are
compositional, so ANY drop in the neuronal denominator inflates every
non-neuronal fraction whether or not a single glial cell was gained.

The fix is to ask a question the confound cannot reach: within the non-neuronal
compartment alone, does the mix change with age? That conditions the neuronal
denominator away entirely.

Result: the microglia / pericyte / smooth-muscle increases DO NOT survive.
They were the neuronal denominator moving. What survives is endothelium and
border-associated macrophages DECLINING as a share of non-neuronal cells.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aging_composition import bh_fdr, cliffs_delta, clr, load  # noqa: E402

LABEL = {
    "333 Endo NN": "Endothelial", "331 Peri NN": "Pericyte",
    "332 SMC NN": "Smooth muscle", "330 VLMC NN": "VLMC",
    "334 Microglia NN": "Microglia", "335 BAM NN": "Border macrophage",
    "338 Lymphoid NN": "Lymphoid", "337 DC NN": "Dendritic cell",
    "326 OPC NN": "OPC", "327 Oligo NN": "Oligodendrocyte",
    "319 Astro-TE NN": "Astrocyte (TE)", "318 Astro-NT NN": "Astrocyte (NT)",
}


def is_non_neuronal(class_name: pd.Series) -> np.ndarray:
    c = class_name.astype(str)
    return (c.str.contains(" NN") | c.str.startswith(("30 ", "31 ", "33 ", "34 "))).to_numpy()


def compose(sub: pd.DataFrame, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, meta = [], []
    for (_, rg, ag), s in sub.groupby(
            ["donor_label", "anatomical_division_label", "donor_age_category"],
            observed=True):
        if len(s) < min_cells:
            continue
        rows.append(s.subclass_name.astype(str).value_counts(normalize=True))
        meta.append({"region": rg, "age": ag})
    return pd.DataFrame(rows).fillna(0.0), pd.DataFrame(meta)


def test(P: pd.DataFrame, M: pd.DataFrame) -> pd.DataFrame:
    C = pd.DataFrame(clr(P.to_numpy()), columns=P.columns)
    out = []
    for t in P.columns:
        A, B, nreg = [], [], 0
        for r in M.region.unique():
            m = (M.region == r).to_numpy()
            a, b = m & (M.age == "adult").to_numpy(), m & (M.age == "aged").to_numpy()
            if a.sum() < 3 or b.sum() < 3:
                continue
            if (P.loc[a, t] > 0).mean() < 0.6 or (P.loc[b, t] > 0).mean() < 0.6:
                continue
            base = np.median(C.loc[a, t])
            A.append(C.loc[a, t].to_numpy() - base)
            B.append(C.loc[b, t].to_numpy() - base)
            nreg += 1
        if nreg < 2:
            continue
        A, B = np.concatenate(A), np.concatenate(B)
        out.append({"subclass": t, "label": LABEL.get(t, t), "nreg": nreg,
                    "clr": float(np.median(B) - np.median(A)),
                    "delta": cliffs_delta(B, A),
                    "pct_adult": 100 * np.median(P.loc[(M.age == "adult").to_numpy(), t]),
                    "pct_aged": 100 * np.median(P.loc[(M.age == "aged").to_numpy(), t]),
                    "p": stats.mannwhitneyu(B, A, alternative="two-sided")[1]})
    R = pd.DataFrame(out)
    R["q"] = bh_fdr(R.p.to_numpy())
    return R.sort_values("clr", ascending=False)


def main() -> int:
    csv = sys.argv[1] if len(sys.argv) > 1 else "aging_cells.csv"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    df = load(csv)
    nn_mask = is_non_neuronal(df.class_name)

    print("=== the confound ===")
    print("  neuronal share of unbiased cells: adult %.1f%%  aged %.1f%%"
          % (100 * (~nn_mask[(df.donor_age_category == "adult").to_numpy()]).mean(),
             100 * (~nn_mask[(df.donor_age_category == "aged").to_numpy()]).mean()))

    whole = test(*compose(df, 150))
    nn = test(*compose(df[nn_mask], 100))

    print("\n=== A. share of ALL cells (confounded by neuronal dropout) ===")
    print("  %-20s %+7s %8s %s" % ("cell type", "CLR", "q", ""))
    for _, r in whole.iterrows():
        if r.subclass not in LABEL:
            continue
        print("  %-20s %+7.3f %8.1e %s"
              % (r.label, r.clr, r.q, "*" if r.q < 0.05 else ""))

    print("\n=== B. share of NON-NEURONAL cells only (confound removed) ===")
    print("  %-20s %+7s %8s %8s %8s %s"
          % ("cell type", "CLR", "%adult", "%aged", "q", ""))
    for _, r in nn.iterrows():
        print("  %-20s %+7.3f %7.2f%% %7.2f%% %8.1e %s"
              % (r.label, r.clr, r.pct_adult, r.pct_aged, r.q,
                 "*" if r.q < 0.05 else ""))

    merged = whole.merge(nn, on="subclass", suffixes=("_whole", "_nn"))
    merged.to_csv(os.path.join(outdir, "results_aging_confound.csv"), index=False)
    nn.to_csv(os.path.join(outdir, "results_aging_nonneuronal.csv"), index=False)

    print("\n=== what changed ===")
    for _, r in merged.iterrows():
        flip = (r.q_whole < 0.05) != (r.q_nn < 0.05)
        if flip:
            print("  %-20s q %.3f -> %.3f   %s"
                  % (r.label_whole, r.q_whole, r.q_nn,
                     "LOST significance" if r.q_whole < 0.05 else "GAINED significance"))
    print("\nwrote results_aging_confound.csv, results_aging_nonneuronal.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
