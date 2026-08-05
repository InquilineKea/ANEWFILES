# Age-heterogeneous human endothelial cell atlas → h5ad

Tooling to pull a human endothelial-cell montage with maximum **age
heterogeneity** as `.h5ad`, ready for downstream statistics.

**The download did not run here.** This sandbox has no egress to any data host
(see [Why nothing was downloaded](#why-nothing-was-downloaded)). What is in
this directory is the tested machinery to run it on a machine that does have
network, plus a synthetic file with the identical schema so analysis code can
be written immediately.

## Files

| File | What it is |
|---|---|
| `fetch_endothelial_atlas.py` | Main script. Queries CZ CELLxGENE Census, writes `.h5ad`. |
| `hsapdv_age.py` | Parses `development_stage` labels → one numeric age axis. Self-testing. |
| `make_demo_h5ad.py` | Generates the synthetic schema-demo file. |
| `endothelial_demo.h5ad` | 1800 cells × 30 genes, 15 donors, ages −0.69 … 89 yr. **Synthetic.** |

## Quick start

```bash
pip install cellxgene-census anndata pandas numpy scipy

python3 fetch_endothelial_atlas.py --self-test                    # offline, ~1s
python3 fetch_endothelial_atlas.py --list-cell-types --tissue lung # no download
python3 fetch_endothelial_atlas.py --tissue lung --max-cells 50000 # small real pull
```

Start with `--tissue lung --max-cells 50000`. An unfiltered all-tissue
endothelial query plausibly returns 1–2M cells and several GB; `--max-cells`
exists so you find that out cheaply rather than mid-download.

## Why CZ CELLxGENE Census

| Source | Age heterogeneity | Verdict |
|---|---|---|
| **CELLxGENE Census** | fetal → 80s+ in one query | **chosen** |
| Tabula Sapiens | adult only, ~20s–70s, 24 donors | good fallback, but adult-only |
| Human Protein Atlas | **none** | rejected |

Human Protein Atlas is out on structural grounds, not convenience: its
single-cell layer is a *consensus across donors*, so per-donor age does not
exist in it. Age heterogeneity is not merely low, it is zero by construction.

Tabula Sapiens is the easiest single file and endothelial cells are cleanly
annotated, but it is adult-only. Since you asked for age heterogeneity
specifically, that is the axis to optimise, so Census wins. If Census is being
difficult, `--local TabulaSapiens.h5ad` runs the same filtering and age
annotation over a file you already have.

## Two traps this handles

### 1. Ontology descendants are not expanded for you

Filtering on `cell_type == 'endothelial cell'` (or on `CL:0000115`) silently
drops capillary, arterial, venous and lymphatic endothelium. Census does **not**
walk the ontology. The script enumerates the labels actually present first
(`--list-cell-types`), then filters on the explicit list.

Corneal endothelium is excluded by default: it matches the string but is
neural-crest derived, not vascular. Override with `--exclude`.

### 2. `development_stage` is a label, not a number

The same quantity appears in ~8 incompatible surface forms. All are parsed to
one monotonic float, **in years, with birth at 0.0 and prenatal ages negative**
(20 pcw → −0.385 yr) so a model spanning gestation to adulthood does not hit a
discontinuity at birth.

Every parse carries a precision tag:

| obs column | meaning |
|---|---|
| `age_years` | the number. `NaN` if unparseable — never silently 0.0 |
| `age_precision` | `exact_year`, `range_midpoint`, `decade_midpoint`, `pcw`, `carnegie_APPROX`, `coarse_*` |
| `age_span_years` | width of the interval the number summarises |
| `age_is_imprecise` | `True` for labels too coarse to use as a point estimate |

`"human adult stage"` parses to 54.5 but is tagged
`coarse_adult_VERY_IMPRECISE`, because it spans ~40 years. Treating that as a
point estimate is the main way this analysis goes wrong quietly, so filter on
`~age_is_imprecise` for anything age-continuous. Unparsed labels are reported
at the end of a run, with counts, so the parser can be extended rather than
guessed at.

## Verification actually performed

Everything below was executed, not assumed:

- **Age parser: 22/22 cases pass**, plus an explicit monotonicity check across
  birth (Carnegie → pcw → newborn → months → years strictly increasing).
- **A real bug was caught by the self-test.** `"45-49 year-old"` contains the
  substring `"49 year-old"`, so the single-year pattern matched the tail of the
  range and returned **49 instead of the midpoint 47**. Ranges are now matched
  first, and the case is pinned in `_CASES`.
- **API signatures verified against the installed library** (`cellxgene-census`
  1.18.0), not recalled. This matters: `get_anndata` takes
  `obs_value_filter` / `obs_column_names`, while `get_obs` takes
  `value_filter` / `column_names`. They are different keywords and mixing them
  raises `TypeError`.
- **`--local` path run end-to-end** on the demo file: 1800 cells, 7 endothelial
  labels, 0 unparsed ages, 9 precision tags.
- **Corneal exclusion fires**; generated filter string checked.
- **Spatial demo is genuinely clustered**: mean nearest-neighbour distance is
  0.54× the CSR expectation, so point-pattern code has real structure to detect.
- **Embedded age signal recoverable**: `corr(age, CDKN2A) = +0.33`,
  `corr(age, LMNB1) = −0.28`, `corr(age, PECAM1) = −0.01`.

Honest caveat on subsampling: `stratified_subsample` provably retains all 15
donors at every budget tested, but on a file this small *naive* uniform
sampling happens to retain them too — so the demo does **not** prove the
stratification earns its keep. It matters at Census scale, where one adult
dataset can outnumber a fetal one 100:1 and uniform sampling flattens the very
axis you came for.

## On "spatial statistics"

Worth deciding before you download. Single-cell Census data is **dissociated**:
there are no tissue coordinates, and `obsm['spatial']` is deliberately left
absent on the real object (with a note in `.uns['spatial_note']`). "Spatial
statistics" on that object means expression-space or organ-level structure
until you join real geometry.

`endothelial_demo.h5ad` *does* ship a synthetic clustered `obsm['spatial']`,
purely so Ripley's K / Moran's I implementations have a non-random pattern to
run against — a CSR null cannot distinguish a correct implementation from one
that always returns zero.

If you want genuine point patterns *with* age variation, that is a different
and much narrower dataset — HuBMAP CODEX / Xenium-style imaging, where donor
age range is far smaller than Census. There is a real trade-off between age
heterogeneity and true spatial coordinates, and no current human resource
gives you both at once. Which one is load-bearing for your analysis?

## Confidence

- **High**: the ontology-descendant problem; the age-label surface forms; the
  HPA and Tabula Sapiens characterisations; all API signatures (checked against
  the installed package).
- **Medium**: exact Census column availability for a given `--tissue` — the
  schema shifts between Census releases; `--list-cell-types` is the cheap probe.
- **Lower**: the Carnegie-stage → post-conception-week table is approximate
  (staging is morphological, sources disagree by a few days) and only matters
  if you keep embryonic samples. Tagged `carnegie_APPROX` for that reason.
- **Unknown until run**: how many cells the unfiltered endothelial query
  returns, and total download size.

## Why nothing was downloaded

Outbound HTTPS from this sandbox goes through a policy-enforcing proxy that
returned **HTTP 403 on `CONNECT`** for every data host tried:

```
api.cellxgene.cziscience.com      403
datasets.cellxgene.cziscience.com 403
www.proteinatlas.org              403
github.com                        400
```

403 on CONNECT is an organisation egress-policy denial — the destination is not
on the allowlist. The proxy documentation is explicit that these must be
reported rather than retried or routed around, so no workaround was attempted.
`WebFetch` is subject to the same policy and also returned 403.

PyPI *is* reachable (it is in `no_proxy`), which is why `cellxgene-census`
could be installed and its API introspected — that is how the signatures above
were verified rather than guessed. To actually download data, run these scripts
somewhere with normal network access, or ask an admin to allowlist
`*.cziscience.com`.
