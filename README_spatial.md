# Real spatial data + cell-type frequencies (Allen Brain Cell Atlas)

This is the follow-up to `README_endothelial.md`, and it supersedes that file's
central claim that **nothing could be downloaded**. Something could.

## What changed

The earlier conclusion — no data reachable — was correct about the hosts tried
and wrong as a generalisation. Re-probing more hosts found that
`s3.amazonaws.com` is **not** blocked, while the portal front-ends are:

| Host | Result |
|---|---|
| `cf.10xgenomics.com`, `portal.hubmapconsortium.org`, `figshare.com` | 403 |
| `census.cziscience.com` (Census release directory) | 403 |
| `raw.githubusercontent.com` | **200** |
| `allen-brain-cell-atlas.s3.us-west-2.amazonaws.com` | **200** |
| `cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com` | **200** |

So the block is on portal/CDN hosts, not on object storage. The Allen Brain
Cell Atlas publishes to an anonymously-readable S3 bucket, and **1.6 GB of real
data was downloaded and analysed.** No credentials, no boto3 — plain HTTPS.

## Data

`MERFISH-C57BL6J-638850-CCF` — whole mouse brain MERFISH:
**3,739,961 cells**, 53 sections, 34 classes / 338 subclasses, with CCF
anatomical parcellation *and* per-cell x/y/z coordinates in millimetres.

```bash
python3 abc_fetch.py --list
python3 abc_fetch.py --dataset merfish      # 1.6 GB
python3 abc_frequencies.py merfish_cells.csv
```

Outputs (committed here):

| File | Contents |
|---|---|
| `abc_hippocampus.h5ad` | 304,642 HPF cells, real `obsm['spatial']` + `spatial_3d` |
| `abc_endothelial.h5ad` | 263,210 endothelial cells brain-wide, real coordinates |
| `results_hpf_subclass_pct.csv` | subclass % per hippocampal subregion |
| `results_endothelial_by_division.csv` | endothelial share per brain division |

## Cell-type frequencies in hippocampus

304,642 cells in the hippocampal formation. Top subclasses:

| Subclass | Share | n |
|---|---|---|
| 319 Astro-TE NN | 19.83% | 60,420 |
| 037 DG Glut | 11.11% | 33,839 |
| 016 CA1-ProS Glut | 7.83% | 23,842 |
| **333 Endo NN** | **7.29%** | **22,206** |
| 327 Oligo NN | 5.79% | 17,637 |

Dominant type per subregion — the anatomy comes out correctly, which is the
best available check that the labels and parcellation are joined properly:

| Region | n | Dominant subclass | Share |
|---|---|---|---|
| DG | 51,599 | 037 DG Glut | 45.85% |
| CA3 | 47,013 | 017 CA3 Glut | 27.53% |
| CA1 | 56,402 | 016 CA1-ProS Glut | 38.51% |
| CA2 | 3,173 | 319 Astro-TE NN | 25.40% |
| ENTm | 39,486 | 008 L2/3 IT ENT Glut | 22.13% |

CA2 is the one that looks wrong at a glance and is not: CA2 is a thin
transitional field with only 3,173 segmented cells, its pyramidal band is
narrow, so astrocytes out-number the projection neurons within the parcellation
boundary. Treat small-region frequencies as boundary-sensitive.

## Endothelial frequencies

Brain-wide vascular composition:

| Subclass | Type | n | % of all cells |
|---|---|---|---|
| 333 Endo NN | endothelial | 263,210 | 7.04% |
| 330 VLMC NN | leptomeningeal | 70,506 | 1.89% |
| 331 Peri NN | pericyte | 46,134 | 1.23% |
| 332 SMC NN | smooth muscle | 26,699 | 0.71% |
| 329 ABC NN | arachnoid barrier | 19,081 | 0.51% |

Endothelial share varies ~2× across brain divisions — hindbrain highest,
cortex/hippocampus middling:

| Division | Endothelial share |
|---|---|
| MY (medulla) | 11.29% |
| P (pons) | 10.05% |
| MB (midbrain) | 8.57% |
| Isocortex | 7.89% |
| HPF | 7.29% |
| HY (hypothalamus) | 6.46% |

Within HPF, endothelial share runs 5.1% (SUB) to 12.4% (IG). **Do not read this
as vascular density.** It is a *proportion of segmented cells*, so it moves when
the neuronal denominator moves. DG looks vessel-poor (5.79%) largely because the
granule layer is extraordinarily cell-dense. A density claim needs vessels per
unit volume, which requires the coordinates, not the frequency table.

## A first spatial result

The point of getting real coordinates was to do spatial statistics. One ran
already, comparing mean nearest-neighbour distance against a random subset of
identical size drawn from the same section (ratio < 1 = clustered, ≈ 1 = random):

| Subclass | n | NN ratio | Interpretation |
|---|---|---|---|
| 037 DG Glut | 5,001 | **0.33** | strongly clustered — the granule blade |
| 016 CA1-ProS Glut | 2,559 | **0.45** | clustered — pyramidal layer |
| 333 Endo NN | 1,897 | **1.05** | indistinguishable from random |

This is a genuine result, not just a smoke test. Neurons are laminated;
endothelium is spatially near-uniform, which is what vasculature must do to
perfuse tissue evenly. It also validates the coordinates: laminar structure
this sharp cannot arise from scrambled positions.

Caveat: this compares against a random subsample of *observed cells*, so it
controls for tissue outline but not for local density. A proper Ripley's K with
an edge-corrected, tissue-shaped null is the next step.

## The age/space trade-off, now settled empirically

Last round this was a judgement call. It is now measured:

| Dataset | Coordinates | Donors | Age range |
|---|---|---|---|
| MERFISH C57BL6J-638850 | **yes** | **1** | none |
| Xenium marmoset (HMBA-CJ23.56.004-BG) | **yes** | **1** | none (5.7 yr, single male) |
| Zeng-Aging-Mouse-10Xv3 | no | **108** | 2 groups: ~P56 vs ~P550 |

**You cannot get both from these resources.** The Xenium dataset you asked
about is real and downloadable — but it is one marmoset, basal ganglia, no
hippocampus, and age heterogeneity is exactly zero, so it cannot support any
age question. The aging dataset has 108 donors and a genuine young/aged
contrast, and no coordinates at all.

If age is the real variable, `aging-mouse` is the dataset and spatial statistics
are off the table. If spatial organisation is the real variable, MERFISH is
excellent and age is off the table. Worth choosing deliberately.

## Statistical caveats that actually bite

1. **Single donor.** Every frequency here describes one mouse. The 53 sections
   are pseudoreplicates: they cannot support a confidence interval that
   generalises to mice. `section_dispersion()` reports median/IQR across
   sections and is labelled dispersion, not standard error, for this reason.
   The spread is not small — CA1's dominant subclass ranges 27.8%–43.3% across
   sections purely as a function of where you cut.

2. **Frequencies are compositional.** They sum to 1 within a region, so they are
   not independent — one type rising forces others down. Ordinary correlations
   or t-tests on raw proportions manufacture spurious negative associations.
   Use `clr_transform()` (provided) before any such test.

3. **Proportion ≠ density.** As above. The denominator is segmented cells, and
   it varies hugely by region.

4. **Segmentation-dependent.** MERFISH cell counts depend on the segmentation;
   small cells in dense layers are under-recovered. Cross-region comparisons
   inherit that bias.

5. **Mouse, not human.** This is the resource that was reachable and that has
   hippocampal parcellation. Human vascular fractions differ.
