#!/usr/bin/env python3
"""
Generate a small SYNTHETIC endothelial h5ad with the same schema the real
Census pull produces, so downstream analysis code can be written and tested
before a multi-GB download finishes.

    python3 make_demo_h5ad.py -o endothelial_demo.h5ad

The expression values are noise. Do not draw biological conclusions from this
file -- it exists to pin the schema and to give spatial-statistics code a
non-CSR point pattern to run against.

Unlike the real single-cell object, this demo DOES populate obsm['spatial']
with deliberately clustered synthetic coordinates, so Ripley's K / Moran's I
implementations have something with real structure to detect (a complete
spatial randomness null would make a broken implementation look correct).
"""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from hsapdv_age import parse_development_stage

# Deliberately spans the full range and exercises every parser branch,
# including the coarse labels you will want to exclude.
DONORS = [
    ("D01", "Carnegie stage 13", "male"),
    ("D02", "8 post-conception weeks human stage", "female"),
    ("D03", "20 post-conception week human stage", "male"),
    ("D04", "newborn human stage", "female"),
    ("D05", "6-month-old human stage", "male"),
    ("D06", "2-year-old human stage", "female"),
    ("D07", "adolescent stage", "male"),
    ("D08", "third decade human stage", "female"),
    ("D09", "34-year-old human stage", "male"),
    ("D10", "45-49 year-old stage", "female"),
    ("D11", "57-year-old human stage", "male"),
    ("D12", "65-69 year-old stage", "female"),
    ("D13", "eighth decade human stage", "male"),
    ("D14", "89-year-old human stage", "female"),
    ("D15", "human adult stage", "male"),  # coarse: should be flagged
]

CELL_TYPES = [
    "endothelial cell of vascular tree",
    "capillary endothelial cell",
    "endothelial cell of artery",
    "vein endothelial cell",
    "endothelial cell of lymphatic vessel",
    "lung microvascular endothelial cell",
    "endothelial cell of sinusoid",
]

GENES = [
    "PECAM1", "CDH5", "VWF", "CLDN5", "ERG", "FLT1", "KDR", "TEK", "ENG",
    "PROX1", "PDPN", "LYVE1", "GJA5", "HEY2", "NR2F2", "ACKR1", "CA4",
    "APLNR", "ESM1", "CXCR4", "SELE", "ICAM1", "VCAM1", "NOS3", "THBD",
    "CDKN1A", "CDKN2A", "GLB1", "LMNB1", "TP53",
]


def clustered_points(n: int, rng: np.random.Generator, n_centers: int = 12) -> np.ndarray:
    """Thomas-process-ish clustered point pattern in a unit square.

    Clustered on purpose: a CSR (complete spatial randomness) pattern cannot
    distinguish a correct Ripley's K from one that always returns ~0.
    """
    centers = rng.uniform(0.05, 0.95, size=(n_centers, 2))
    which = rng.integers(0, n_centers, size=n)
    pts = centers[which] + rng.normal(scale=0.035, size=(n, 2))
    return np.clip(pts, 0.0, 1.0)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output", default="endothelial_demo.h5ad")
    p.add_argument("--cells-per-donor", type=int, default=120)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    n_per = args.cells_per_donor
    n = n_per * len(DONORS)
    g = len(GENES)

    donor_id, stage, sex = [], [], []
    for d, s, sx in DONORS:
        donor_id += [d] * n_per
        stage += [s] * n_per
        sex += [sx] * n_per

    obs = pd.DataFrame(
        {
            "soma_joinid": np.arange(n),
            "dataset_id": "SYNTHETIC-DEMO",
            "assay": "10x 3' v3",
            "cell_type": rng.choice(CELL_TYPES, size=n),
            "cell_type_ontology_term_id": "CL:0000115",
            "development_stage": stage,
            "disease": "normal",
            "donor_id": donor_id,
            "sex": sex,
            "suspension_type": "cell",
            "tissue": "lung",
            "tissue_general": "lung",
            "is_primary_data": True,
        },
        index=pd.Index(["cell_%05d" % i for i in range(n)], name="cell_id"),
    )

    parsed = [parse_development_stage(x) for x in obs["development_stage"]]
    obs["age_years"] = np.array([q.years for q in parsed], dtype=float)
    obs["age_precision"] = pd.Categorical([q.precision for q in parsed])
    obs["age_span_years"] = np.array([q.span_years for q in parsed], dtype=float)
    obs["age_is_imprecise"] = np.array([q.is_imprecise for q in parsed], dtype=bool)
    for c in ("cell_type", "donor_id", "sex", "development_stage", "tissue", "tissue_general"):
        obs[c] = obs[c].astype("category")

    # Counts with a weak, genuinely present age signal on the senescence-ish
    # genes, so a correct age regression finds something and a broken one does not.
    age = obs["age_years"].to_numpy()
    age_z = np.nan_to_num((age - np.nanmean(age)) / (np.nanstd(age) + 1e-9))
    base = rng.gamma(shape=2.0, scale=1.5, size=(n, g))
    for gene, direction in (("CDKN1A", 1.0), ("CDKN2A", 1.0), ("LMNB1", -1.0)):
        j = GENES.index(gene)
        base[:, j] *= np.exp(0.45 * direction * age_z)
    counts = rng.poisson(base).astype(np.float32)
    counts[rng.random(counts.shape) < 0.35] = 0.0  # dropout

    var = pd.DataFrame(
        {"feature_name": GENES, "feature_id": ["ENSGDEMO%05d" % i for i in range(g)]},
        index=pd.Index(GENES, name="gene"),
    )

    adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.obsm["spatial"] = clustered_points(n, rng)
    adata.uns["endothelial_query"] = {
        "source": "SYNTHETIC DEMO -- not real data",
        "age_convention": "years relative to term birth (40 pcw); prenatal negative",
    }
    adata.uns["spatial_note"] = (
        "obsm['spatial'] here is SYNTHETIC and clustered, for testing point-pattern "
        "code only. The real Census pull has no coordinates and omits this slot."
    )

    adata.write_h5ad(args.output, compression="gzip")
    print("wrote %s: %d cells x %d genes, %d donors"
          % (args.output, adata.n_obs, adata.n_vars, obs["donor_id"].nunique()))
    ok = obs["age_years"].notna()
    print("age range: %.3f .. %.3f years (%d/%d parsed, %d imprecise)"
          % (obs.loc[ok, "age_years"].min(), obs.loc[ok, "age_years"].max(),
             int(ok.sum()), n, int(obs["age_is_imprecise"].sum())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
