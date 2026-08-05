#!/usr/bin/env python3
"""
FTU innervation, evaluated on Allen Brain Cell Atlas MERFISH.

    python3 abc_fetch.py --dataset merfish
    python3 ftu_innervation.py merfish_cells.csv

Scope, and one thing this CANNOT do
-----------------------------------
"FTU" (functional tissue unit) is HuBMAP vocabulary — glomerulus, alveolus,
crypt. No HuBMAP or Human Reference Atlas host is reachable from this
environment (portal, assets and cdn.humanatlas.io all fail), so human
peripheral-organ FTUs are unavailable. What is available is the brain, where
the canonical functional units are the hippocampal subfields — DG, CA3, CA2,
CA1 and the subicular/entorhinal outputs form the trisynaptic circuit.

**Soma-based spatial transcriptomics cannot measure innervation.** Innervation
is axons and terminals. MERFISH segments and types cell BODIES. A cholinergic
axon arborising densely through CA1 leaves no cholinergic cell body in CA1 —
its soma is in the medial septum, millimetres away. So counting neuromodulatory
somata inside an FTU measures almost nothing about how strongly that FTU is
innervated, and any script that reports "cholinergic cells per FTU" as an
innervation score is measuring the wrong object.

This script therefore reports three things it CAN establish:

  1. **Intrinsic circuit composition** per FTU — excitatory/inhibitory balance
     and interneuron diversity. This is local and genuinely measurable.
  2. **That modulatory innervation is extrinsic**, quantitatively: how many
     modulatory somata lie inside each FTU (≈0) and how far the nearest one is.
  3. **The source nuclei** that must therefore supply it, with their sizes —
     the upstream half of the circuit, which soma data does capture well.

The projection itself — which source reaches which FTU, and how densely — is
exactly what this modality cannot see, and needs tract tracing or axonal
labelling.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

MODULATORS = ["Chol", "Dopa", "Sero", "Nora", "Hist"]
FAST_NT = ["Glut", "GABA", "GABA-Glyc", "Glut-GABA"]

# Hippocampal formation subfields, in trisynaptic order. These are the brain's
# canonical functional tissue units.
FTUS = ["DG", "CA3", "CA2", "CA1", "SUB", "ProS", "POST", "PRE", "PAR",
        "ENTl", "ENTm"]

COLS = ["neurotransmitter", "class", "subclass", "parcellation_division",
        "parcellation_structure", "x_reconstructed", "y_reconstructed",
        "z_reconstructed"]


def shannon(counts: np.ndarray) -> float:
    c = np.asarray(counts, float)
    c = c[c > 0]
    if c.sum() <= 0:
        return np.nan
    p = c / c.sum()
    return float(-(p * np.log(p)).sum())


def load(csv: str) -> pd.DataFrame:
    cache = os.path.splitext(csv)[0] + "_ftu.pkl"
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    df = pd.read_csv(csv, usecols=COLS,
                     dtype={c: "category" for c in COLS if not c.startswith(("x_", "y_", "z_"))})
    df.to_pickle(cache)
    return df


def intrinsic_profile(df: pd.DataFrame) -> pd.DataFrame:
    """Excitatory/inhibitory balance and interneuron diversity per FTU.

    Non-neuronal cells carry no neurotransmitter label and are excluded from
    the E/I denominator -- otherwise the ratio would track glial abundance
    rather than circuit composition.
    """
    nt = df.neurotransmitter.astype(str)
    st = df.parcellation_structure.astype(str)
    neurons = nt.isin(FAST_NT) | nt.isin(MODULATORS)
    rows = []
    for f in FTUS:
        m = (st == f).to_numpy()
        if m.sum() < 200:
            continue
        sub_nt = nt[m & neurons.to_numpy()]
        n_neu = len(sub_nt)
        if n_neu < 100:
            continue
        glut = int((sub_nt == "Glut").sum())
        gaba = int(sub_nt.isin(["GABA", "GABA-Glyc"]).sum())
        # diversity of GABAergic subclasses = interneuron repertoire
        gsub = df.subclass.astype(str)[m & nt.isin(["GABA", "GABA-Glyc"]).to_numpy()]
        _, gc = np.unique(gsub, return_counts=True)
        rows.append({
            "FTU": f, "n_cells": int(m.sum()), "n_neurons": n_neu,
            "pct_glut": 100 * glut / n_neu, "pct_gaba": 100 * gaba / n_neu,
            "EI_ratio": glut / max(gaba, 1),
            "n_gaba_subclasses": len(gc),
            "interneuron_entropy": shannon(gc),
            "n_modulatory_somata": int(sub_nt.isin(MODULATORS).sum()),
        })
    return pd.DataFrame(rows)


def extrinsic_distance(df: pd.DataFrame) -> pd.DataFrame:
    """Distance from each FTU's cells to the nearest modulatory soma.

    Uses 3D reconstructed coordinates (mm). If innervation were locally
    sourced these distances would be tens of microns; if extrinsic they are
    millimetres.
    """
    nt = df.neurotransmitter.astype(str)
    st = df.parcellation_structure.astype(str)
    xyz = df[["x_reconstructed", "y_reconstructed", "z_reconstructed"]].to_numpy(float)
    ok = np.isfinite(xyz).all(axis=1)

    trees = {}
    for mod in MODULATORS:
        m = (nt == mod).to_numpy() & ok
        if m.sum() >= 10:
            trees[mod] = cKDTree(xyz[m])

    rows = []
    for f in FTUS:
        m = (st == f).to_numpy() & ok
        if m.sum() < 200:
            continue
        pts = xyz[m]
        rec = {"FTU": f}
        for mod, tree in trees.items():
            d, _ = tree.query(pts, k=1)
            rec["d_%s_median_mm" % mod] = float(np.median(d))
            rec["d_%s_min_mm" % mod] = float(np.min(d))
        rows.append(rec)
    return pd.DataFrame(rows)


def source_nuclei(df: pd.DataFrame, top: int = 5) -> pd.DataFrame:
    """Where the modulatory somata actually are -- the upstream half."""
    nt = df.neurotransmitter.astype(str)
    st = df.parcellation_structure.astype(str)
    dv = df.parcellation_division.astype(str)
    rows = []
    for mod in MODULATORS:
        m = (nt == mod).to_numpy()
        if m.sum() == 0:
            continue
        vc = st[m].value_counts().head(top)
        for struct, n in vc.items():
            rows.append({"modulator": mod, "n_total": int(m.sum()),
                         "source_structure": struct,
                         "division": dv[m][st[m] == struct].mode().iat[0],
                         "n_somata": int(n),
                         "pct_of_modulator": 100 * n / m.sum()})
    return pd.DataFrame(rows)


def main() -> int:
    csv = sys.argv[1] if len(sys.argv) > 1 else "merfish_cells.csv"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    if not os.path.exists(csv):
        print("missing %s -- run: python3 abc_fetch.py --dataset merfish" % csv)
        return 2
    df = load(csv)
    print("loaded %d cells" % len(df))

    print("\n" + "=" * 78)
    print("1. INTRINSIC CIRCUIT COMPOSITION PER FTU  (measurable from somata)")
    print("=" * 78)
    I = intrinsic_profile(df)
    print("  %-6s %9s %9s %8s %8s %7s %8s %6s"
          % ("FTU", "cells", "neurons", "%Glut", "%GABA", "E/I", "IN-types", "H"))
    for _, r in I.iterrows():
        print("  %-6s %9d %9d %7.1f%% %7.1f%% %7.1f %8d %6.2f"
              % (r.FTU, r.n_cells, r.n_neurons, r.pct_glut, r.pct_gaba,
                 r.EI_ratio, r.n_gaba_subclasses, r.interneuron_entropy))

    print("\n" + "=" * 78)
    print("2. IS MODULATORY INNERVATION LOCAL?  (modulatory somata inside FTU)")
    print("=" * 78)
    D = extrinsic_distance(df)
    M = I[["FTU", "n_neurons", "n_modulatory_somata"]].merge(D, on="FTU")
    dcols = [c for c in M.columns if c.startswith("d_") and c.endswith("median_mm")]
    print("  %-6s %10s %10s   %s"
          % ("FTU", "mod.somata", "per 10k", "median distance to nearest soma (mm)"))
    print("  %-6s %10s %10s   %s"
          % ("", "", "", "  ".join("%-6s" % c.split("_")[1] for c in dcols)))
    for _, r in M.iterrows():
        per10k = 1e4 * r.n_modulatory_somata / r.n_neurons
        print("  %-6s %10d %10.1f   %s"
              % (r.FTU, r.n_modulatory_somata, per10k,
                 "  ".join("%-6.2f" % r[c] for c in dcols)))

    print("\n" + "=" * 78)
    print("3. SOURCE NUCLEI  (where the innervation must originate)")
    print("=" * 78)
    S = source_nuclei(df)
    for mod in MODULATORS:
        s = S[S.modulator == mod]
        if not len(s):
            continue
        print("\n  %-5s  %d somata brain-wide" % (mod, s.n_total.iat[0]))
        for _, r in s.iterrows():
            print("      %-18s %-8s %6d  (%4.1f%%)"
                  % (r.source_structure, r.division, r.n_somata, r.pct_of_modulator))

    I.to_csv(os.path.join(outdir, "results_ftu_intrinsic.csv"), index=False)
    M.to_csv(os.path.join(outdir, "results_ftu_extrinsic.csv"), index=False)
    S.to_csv(os.path.join(outdir, "results_ftu_sources.csv"), index=False)
    print("\nwrote results_ftu_{intrinsic,extrinsic,sources}.csv")
    print("\nNOTE: distances above are soma-to-soma. They establish that modulatory")
    print("input is extrinsic; they do NOT quantify innervation density, which")
    print("requires axonal labelling this modality does not provide.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
