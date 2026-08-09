#!/usr/bin/env python3
"""
External check: my age-composition results vs Allen's own published age analysis.

`cluster.csv` in the Zeng-Aging-Mouse release carries the data producers' own
per-cluster adult/aged counts and odds ratios. Aggregating those to subclass
gives an independent yardstick computed by different people with a different
method (raw odds ratios, no CLR, no region stratification).

    python3 validate_against_allen.py aging_cluster.csv

Result: Spearman rho = +0.80 against my all-cells analysis and +0.78 against
the non-neuronal one, so the pipeline reproduces their ranking. Allen's own
numbers also put endothelium (-0.200) and border macrophages (-0.410) on the
DECLINING side, which matches the corrected finding rather than the retracted
one.

IMPORTANT -- what this does NOT establish. Both analyses run on the same cells
from the same dissociations. A technical artefact in cell recovery is shared,
so agreement here cannot distinguish biology from dissociation bias. This is
concordance of method, not independent replication.
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
from scipy import stats

def main() -> int:
    cl = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "aging_cluster.csv")
    piv = pd.read_csv("cluster_mapping_pivot.csv", usecols=["cluster_alias", "subclass_name"])
    m = cl.merge(piv, on="cluster_alias", how="left")
    g = m.groupby("subclass_name", observed=True)[
        ["number_of_adult_cells", "number_of_aged_cells"]].sum()
    ta, tg = g.number_of_adult_cells.sum(), g.number_of_aged_cells.sum()
    g["allen_log2OR"] = np.log2(((g.number_of_aged_cells + 1) / (tg - g.number_of_aged_cells)) /
                                ((g.number_of_adult_cells + 1) / (ta - g.number_of_adult_cells)))
    mine = pd.read_csv("results_aging_confound.csv")
    c = mine.merge(g.reset_index(), left_on="subclass", right_on="subclass_name")
    for col, name in [("clr_whole", "all-cells"), ("clr_nn", "non-neuronal")]:
        rho, p = stats.spearmanr(c.allen_log2OR, c[col])
        print("Spearman(Allen log2OR, %-13s) = %+.3f  p=%.4f" % (name, rho, p))
    c[["label_whole", "allen_log2OR", "clr_whole", "clr_nn", "q_whole", "q_nn"]] \
        .sort_values("allen_log2OR", ascending=False) \
        .to_csv("results_validation_vs_allen.csv", index=False)
    print("wrote results_validation_vs_allen.csv")
    return 0

if __name__ == "__main__":
    sys.exit(main())
