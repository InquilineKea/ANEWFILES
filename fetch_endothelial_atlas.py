#!/usr/bin/env python3
"""
Fetch a human endothelial-cell montage with maximal age heterogeneity, as .h5ad.

Source: CZ CELLxGENE Census (`cellxgene-census`), which is the only easy human
reference source that spans fetal -> 80s in a single query. See README.md for
why not Human Protein Atlas (no per-donor age at all) and why not Tabula
Sapiens alone (adult-only, ~20s-70s).

Usage
-----
    pip install cellxgene-census anndata pandas numpy scipy

    # always start here -- offline, no network, ~1 second
    python3 fetch_endothelial_atlas.py --self-test

    # see which endothelial labels exist before pulling any expression data
    python3 fetch_endothelial_atlas.py --list-cell-types --tissue lung

    # a small real pull
    python3 fetch_endothelial_atlas.py --tissue lung --max-cells 50000

    # everything, all tissues (large -- expect multiple GB)
    python3 fetch_endothelial_atlas.py --max-cells 500000 -o endo_all.h5ad

Two things this script handles that will otherwise bite you
-----------------------------------------------------------
1. Ontology descendants. Filtering on `cell_type == 'endothelial cell'` (or on
   CL:0000115) silently drops capillary, arterial, venous and lymphatic
   endothelium -- Census does NOT expand the ontology for you. This script
   enumerates the actual labels present first, then filters on the explicit
   list. Corneal endothelium is excluded by default: it matches the string but
   is neural-crest derived, not vascular.

2. Age. `development_stage` is a label, not a number, in ~8 incompatible
   surface forms. See hsapdv_age.py. Everything is parsed to one monotonic
   float with birth at 0.0, tagged by precision, and rows that cannot be
   parsed keep NaN and are REPORTED rather than dropped.
"""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

import numpy as np
import pandas as pd

from hsapdv_age import parse_development_stage, run_self_test

# Labels that contain "endothelial" but are not vascular/lymphatic endothelium.
DEFAULT_EXCLUDE = ("corneal endothelial cell",)

# Fallback list, used only by --local when we cannot enumerate from Census.
KNOWN_ENDOTHELIAL_LABELS = (
    "endothelial cell",
    "endothelial cell of vascular tree",
    "blood vessel endothelial cell",
    "endothelial cell of artery",
    "capillary endothelial cell",
    "vein endothelial cell",
    "endothelial cell of lymphatic vessel",
    "endothelial cell of sinusoid",
    "glomerular endothelial cell",
    "cardiac endothelial cell",
    "lung microvascular endothelial cell",
    "pulmonary artery endothelial cell",
    "endothelial tip cell",
)

# Columns worth carrying. Keeping this tight matters -- Census obs is wide and
# the full set roughly doubles the in-memory footprint of a large pull.
OBS_COLUMNS = [
    "soma_joinid",
    "dataset_id",
    "assay",
    "cell_type",
    "cell_type_ontology_term_id",
    "development_stage",
    "development_stage_ontology_term_id",
    "disease",
    "donor_id",
    "self_reported_ethnicity",
    "sex",
    "suspension_type",
    "tissue",
    "tissue_general",
    "is_primary_data",
]


def _quote_list(values: Sequence[str]) -> str:
    """Render a python list literal for a Census value_filter string."""
    return "[" + ", ".join("'%s'" % v.replace("'", r"\'") for v in values) + "]"


def build_obs_filter(
    cell_types: Sequence[str],
    tissue: str | None,
    healthy_only: bool,
    primary_only: bool,
) -> str:
    clauses = ["cell_type in %s" % _quote_list(cell_types)]
    if tissue:
        clauses.append("tissue_general == '%s'" % tissue)
    if healthy_only:
        clauses.append("disease == 'normal'")
    if primary_only:
        # Census contains the same cells more than once across overlapping
        # collections. Without this you get duplicate donors, which quietly
        # inflates n and breaks any per-donor statistic.
        clauses.append("is_primary_data == True")
    return " and ".join(clauses)


def discover_endothelial_labels(
    census, tissue: str | None, exclude: Sequence[str], healthy_only: bool
) -> list[str]:
    """Enumerate endothelial cell_type labels actually present in Census.

    Uses get_obs (metadata only) rather than get_anndata, so this costs a
    metadata scan and no expression data.
    """
    import cellxgene_census

    clauses = []
    if tissue:
        clauses.append("tissue_general == '%s'" % tissue)
    if healthy_only:
        clauses.append("disease == 'normal'")
    value_filter = " and ".join(clauses) if clauses else None

    # NOTE: get_obs takes `value_filter` / `column_names`, whereas get_anndata
    # takes `obs_value_filter` / `obs_column_names`. They are not the same
    # keyword; mixing them up raises TypeError.
    obs = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter=value_filter,
        column_names=["cell_type"],
    )
    labels = sorted(set(obs["cell_type"].astype(str)))
    endo = [x for x in labels if "endothel" in x.lower()]
    kept = [x for x in endo if x.lower() not in {e.lower() for e in exclude}]
    dropped = [x for x in endo if x not in kept]
    if dropped:
        print("  excluded (matches 'endothelial' but not vascular): %s" % ", ".join(dropped))
    return kept


def annotate_age(obs: pd.DataFrame) -> pd.DataFrame:
    """Add age_years / age_precision / age_span_years columns in place."""
    parsed = [parse_development_stage(x) for x in obs["development_stage"]]
    obs["age_years"] = np.array([p.years for p in parsed], dtype=float)
    obs["age_precision"] = pd.Categorical([p.precision for p in parsed])
    obs["age_span_years"] = np.array([p.span_years for p in parsed], dtype=float)
    obs["age_is_imprecise"] = np.array([p.is_imprecise for p in parsed], dtype=bool)
    return obs


def report_age_coverage(obs: pd.DataFrame) -> None:
    n = len(obs)
    age = obs["age_years"]
    n_bad = int(age.isna().sum())
    print("\n--- age coverage ---")
    print("cells:                %d" % n)
    print("unparsed age:         %d (%.2f%%)" % (n_bad, 100.0 * n_bad / max(n, 1)))
    if n - n_bad:
        print("age range (years):    %.3f .. %.3f" % (age.min(), age.max()))
        print("donors:               %d" % obs["donor_id"].nunique())
        n_imp = int(obs["age_is_imprecise"].sum())
        print("imprecise-age cells:  %d (%.2f%%)  <- consider excluding" % (
            n_imp, 100.0 * n_imp / max(n, 1)))
        print("\nprecision tags:")
        for tag, cnt in obs["age_precision"].value_counts().items():
            print("  %-38s %8d" % (tag, cnt))
    if n_bad:
        bad = obs.loc[age.isna(), "development_stage"].astype(str).value_counts()
        print("\nunparsed development_stage labels (add these to hsapdv_age.py):")
        for lab, cnt in bad.head(20).items():
            print("  %-50r %6d" % (lab, cnt))


def stratified_subsample(
    obs: pd.DataFrame, max_cells: int, seed: int = 0
) -> np.ndarray:
    """Pick <=max_cells row positions, stratified by donor.

    Uniform subsampling would flatten the age axis -- the biggest datasets
    dominate and they are adult-only. Stratifying by donor preserves the age
    spread, which is the whole reason for choosing this source.
    """
    if len(obs) <= max_cells:
        return np.arange(len(obs))
    rng = np.random.default_rng(seed)
    groups = obs.groupby("donor_id", observed=True).indices
    n_groups = len(groups)
    per = max(1, max_cells // max(n_groups, 1))
    picked: list[np.ndarray] = []
    leftover: list[np.ndarray] = []
    for idx in groups.values():
        idx = np.asarray(idx)
        if len(idx) <= per:
            picked.append(idx)
        else:
            sel = rng.choice(idx, size=per, replace=False)
            picked.append(sel)
            leftover.append(np.setdiff1d(idx, sel))
    out = np.concatenate(picked) if picked else np.array([], dtype=int)
    # Top back up to the budget from the donors that had spare cells.
    if len(out) < max_cells and leftover:
        pool = np.concatenate(leftover)
        need = min(max_cells - len(out), len(pool))
        out = np.concatenate([out, rng.choice(pool, size=need, replace=False)])
    return np.sort(out)


def run_census(args: argparse.Namespace) -> int:
    import cellxgene_census

    print("opening Census (version=%s) ..." % args.census_version)
    with cellxgene_census.open_soma(census_version=args.census_version) as census:
        exclude = tuple(args.exclude) if args.exclude else DEFAULT_EXCLUDE

        print("enumerating endothelial cell types ...")
        cell_types = discover_endothelial_labels(
            census, args.tissue, exclude, not args.include_disease
        )
        if not cell_types:
            print("ERROR: no endothelial cell types found for that filter.", file=sys.stderr)
            return 2
        print("  found %d endothelial labels:" % len(cell_types))
        for c in cell_types:
            print("    - %s" % c)

        if args.list_cell_types:
            return 0

        obs_filter = build_obs_filter(
            cell_types, args.tissue, not args.include_disease, not args.include_secondary
        )
        print("\nobs filter:\n  %s" % obs_filter)

        # Pull metadata first so we can size the request and subsample BEFORE
        # materialising any expression matrix.
        print("\nscanning matching cells (metadata only) ...")
        obs_preview = cellxgene_census.get_obs(
            census,
            "homo_sapiens",
            value_filter=obs_filter,
            column_names=["soma_joinid", "donor_id", "development_stage"],
        )
        n_total = len(obs_preview)
        print("  %d cells match." % n_total)
        if n_total == 0:
            print("ERROR: filter matched nothing.", file=sys.stderr)
            return 2

        obs_preview = annotate_age(obs_preview)
        report_age_coverage(obs_preview)

        coords = None
        if n_total > args.max_cells:
            keep = stratified_subsample(obs_preview, args.max_cells, args.seed)
            coords = obs_preview["soma_joinid"].to_numpy()[keep]
            print(
                "\nsubsampling %d -> %d cells, stratified by donor."
                % (n_total, len(coords))
            )

        print("\nfetching expression matrix ... (this is the slow part)")
        adata = cellxgene_census.get_anndata(
            census,
            organism="homo_sapiens",
            measurement_name="RNA",
            X_name=args.x_layer,
            obs_value_filter=obs_filter,
            obs_coords=coords,
            obs_column_names=OBS_COLUMNS,
        )

    adata.obs = annotate_age(adata.obs)
    adata.uns["endothelial_query"] = {
        "source": "CZ CELLxGENE Census",
        "census_version": args.census_version,
        "obs_value_filter": obs_filter,
        "cell_types": list(cell_types),
        "excluded_cell_types": list(exclude),
        "X_layer": args.x_layer,
        "age_convention": "years relative to term birth (40 pcw); prenatal negative",
    }
    return _finish(adata, args)


def run_local(args: argparse.Namespace) -> int:
    """Fallback path: annotate an h5ad you already have on disk."""
    import anndata as ad

    print("reading %s ..." % args.local)
    adata = ad.read_h5ad(args.local)
    print("  %d cells x %d genes" % (adata.n_obs, adata.n_vars))

    col = args.cell_type_col
    if col not in adata.obs.columns:
        print("ERROR: --cell-type-col %r not in obs. Available: %s"
              % (col, list(adata.obs.columns)), file=sys.stderr)
        return 2
    exclude = {e.lower() for e in (args.exclude or DEFAULT_EXCLUDE)}
    labels = adata.obs[col].astype(str)
    mask = labels.str.lower().str.contains("endothel") & ~labels.str.lower().isin(exclude)
    print("  keeping %d endothelial cells (%d labels)"
          % (int(mask.sum()), labels[mask].nunique()))
    adata = adata[mask.to_numpy()].copy()

    stage_col = args.stage_col
    if stage_col not in adata.obs.columns:
        print("WARNING: %r not in obs; age columns will be all-NaN." % stage_col)
        adata.obs["development_stage"] = np.nan
    elif stage_col != "development_stage":
        adata.obs["development_stage"] = adata.obs[stage_col]
    adata.obs = annotate_age(adata.obs)
    report_age_coverage(adata.obs)
    adata.uns["endothelial_query"] = {
        "source": "local file: %s" % args.local,
        "age_convention": "years relative to term birth (40 pcw); prenatal negative",
    }
    return _finish(adata, args)


def _finish(adata, args: argparse.Namespace) -> int:
    # Placeholder for tissue coordinates. Single-cell Census data has NO
    # spatial coordinates -- this slot stays empty unless you join geometry
    # from a spatial assay. See README.md.
    if "spatial" not in adata.obsm:
        adata.uns["spatial_note"] = (
            "No tissue coordinates: this is dissociated single-cell data. "
            "obsm['spatial'] is intentionally absent. Join a spatial assay "
            "(e.g. HuBMAP CODEX) before running point-pattern statistics."
        )
    print("\nwriting %s ..." % args.output)
    adata.write_h5ad(args.output, compression="gzip")
    print("done: %d cells x %d genes -> %s" % (adata.n_obs, adata.n_vars, args.output))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Download an age-heterogeneous human endothelial cell atlas as h5ad.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-o", "--output", default="endothelial_atlas.h5ad")
    p.add_argument("--tissue", default=None,
                   help="tissue_general filter, e.g. lung, heart, kidney. Omit for all.")
    p.add_argument("--max-cells", type=int, default=200_000,
                   help="cap on cells; subsampled stratified by donor (default 200000)")
    p.add_argument("--census-version", default="stable")
    p.add_argument("--x-layer", default="raw", choices=["raw", "normalized"])
    p.add_argument("--include-disease", action="store_true",
                   help="keep non-normal samples (default: disease == 'normal' only)")
    p.add_argument("--include-secondary", action="store_true",
                   help="keep is_primary_data == False (default: primary only; "
                        "leaving this off avoids duplicate donors)")
    p.add_argument("--exclude", nargs="*", default=None,
                   help="cell_type labels to drop (default: %s)" % (DEFAULT_EXCLUDE,))
    p.add_argument("--list-cell-types", action="store_true",
                   help="print endothelial labels and exit without downloading")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--self-test", action="store_true",
                   help="run offline parser tests and exit (no network)")
    p.add_argument("--local", default=None,
                   help="skip Census; annotate this local .h5ad instead")
    p.add_argument("--cell-type-col", default="cell_type", help="(--local only)")
    p.add_argument("--stage-col", default="development_stage", help="(--local only)")
    args = p.parse_args(argv)

    if args.self_test:
        ok, total = run_self_test()
        return 0 if ok == total else 1

    if args.local:
        return run_local(args)

    try:
        import cellxgene_census  # noqa: F401
    except ImportError:
        print("ERROR: cellxgene-census not installed.\n"
              "  pip install cellxgene-census\n"
              "Or annotate a file you already have:\n"
              "  python3 fetch_endothelial_atlas.py --local TabulaSapiens.h5ad",
              file=sys.stderr)
        return 2
    return run_census(args)


if __name__ == "__main__":
    sys.exit(main())
