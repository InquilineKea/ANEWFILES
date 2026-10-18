# ANEWFILES — research project index

A directory of research experiments, interactive demonstrations, and analysis code.

**Across GitHub:** [Browse all 66 public repositories by subject](REPOSITORIES.md) · [Repository inventory JSON](repository-catalog.json)

**Start here:** projects currently live on separate branches. The default branch is an index and archive; cloning it does not put every project into your working directory.

## Projects

| Project | What it contains | Source branch |
|---|---|---|
| [Gramian Lab](https://github.com/InquilineKea/ANEWFILES/blob/claude/gramian-interactive-simulation-s90ud9/gramian-lab/index.html) | Interactive restoration/control Gramian demonstration. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/gramian-interactive-simulation-s90ud9) |
| [Noggin Island](https://github.com/InquilineKea/ANEWFILES/blob/claude/zealous-ptolemy-5ta9rh/noggin-island/index.html) | Browser project. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/zealous-ptolemy-5ta9rh) |
| [Computational matter zoo](https://github.com/InquilineKea/ANEWFILES/blob/claude/gifted-meitner-cppi8d/computational-matter-zoo/REPORT.md) | Report and numerical checks for computational matter. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/gifted-meitner-cppi8d) |
| [Kagome research](https://github.com/InquilineKea/ANEWFILES/blob/claude/nice-cerf-umzz6k/analysis/README.md) | Moiré/Weyl visualizations, claim checks, and AI-Scientist tasks. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/nice-cerf-umzz6k) |
| [Brain and endothelial atlas](https://github.com/InquilineKea/ANEWFILES/blob/claude/endothelial-atlas-download-7vk2je/README_endothelial.md) | Atlas acquisition, brain aging, spatial entropy, and innervation analyses. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/endothelial-atlas-download-7vk2je) |
| [Boston environment](https://github.com/InquilineKea/ANEWFILES/blob/claude/explore-boston-data-mcp-UuLpd/boston_311_pollution.py) | 311 complaints, emissions, demographics, maps, and correlations. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/explore-boston-data-mcp-UuLpd) |
| [GTEx eQTL queries](https://github.com/InquilineKea/ANEWFILES/blob/claude/check-gtex-eqtl-variants-em0bA/check_gtex_eqtl.py) | FOXO3-neighborhood eQTL query script and saved results. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/check-gtex-eqtl-variants-em0bA) |
| [Elegance evaluation](https://github.com/InquilineKea/ANEWFILES/blob/claude/elegant-noether-tt38nx/chatgpt-elegance-eval/README.md) | Conversation evaluation, source transcript, and concept map. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/elegant-noether-tt38nx) |
| [CR Society archive](https://github.com/InquilineKea/ANEWFILES/blob/claude/curl-forum-thread-umkeC/crsociety_forum_thread.txt) | Archived forum text. | [Browse](https://github.com/InquilineKea/ANEWFILES/tree/claude/curl-forum-thread-umkeC) |

## Open or run a project

For an HTML project, follow its link, download the raw HTML, and open the downloaded file in your browser. A GitHub source link is not a deployed app.

For a Python project, make a separate checkout so its files and installed packages stay separate. Example:

```bash
git clone --single-branch --branch claude/nice-cerf-umzz6k https://github.com/InquilineKea/ANEWFILES.git kagome-research
cd kagome-research
python3 -m venv .venv
# macOS/Linux:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r agsci_tasks/kagome_flatband/requirements.txt
python analysis/check_claims.py --quick
```

These are setup instructions derived from the source; this indexing change did not rerun the analyses.

| Project | Setup / entry point |
|---|---|
| Gramian Lab | Browser HTML; download the file and open it. Google Fonts needs internet. |
| Noggin Island | HTML entry point; runtime dependencies have not been audited. |
| Computational matter zoo | Python + NumPy; run `python3 computational-matter-zoo/checks.py`. |
| Kagome research | Claim checks use NumPy; task requirements specify `numpy>=2.0`. See task setup below. |
| Brain and endothelial atlas | See the separate analysis READMEs; the endothelial downloader documents `cellxgene-census anndata pandas numpy scipy`. |
| Boston environment | Python scripts and saved outputs; a complete dependency manifest is still needed. |
| GTEx eQTL queries | Python standard-library imports in the inspected entry point; live queries require network access. |
| Elegance evaluation | Read the Markdown or download `chatgpt-elegance-eval/concept-map.html` and open it. |
| CR Society archive | Plain text. |

The [AI-Scientist task README](https://github.com/InquilineKea/ANEWFILES/blob/claude/nice-cerf-umzz6k/agsci_tasks/README.md) describes its separate AI-Scientist-v3 checkout and Docker setup. Installing the task's NumPy requirement alone does not install the agent framework.

## Scientific status

- **Kagome:** [claim checks](https://github.com/InquilineKea/ANEWFILES/blob/claude/nice-cerf-umzz6k/analysis/README.md) document supported, partial, and failed claims. Read these alongside the visualizations.
- **Brain aging:** [the aging report](https://github.com/InquilineKea/ANEWFILES/blob/claude/endothelial-atlas-download-7vk2je/README_aging.md) records a denominator-confound correction. Earlier conclusions should not be repeated without that correction.
- **Endothelial downloader:** its README identifies `endothelial_demo.h5ad` as synthetic and says that particular download did not run in the original sandbox. Other atlas analyses have their own provenance.
- **Interactive demonstrations:** availability of an HTML file is not evidence of biological or physical validation.
- **Saved results:** reproducibility depends on source data, environment, parameters, and the corresponding code revision.

## Organization status

The project index was added on 2026-10-18. Existing project branches are retained. Separate destination repositories have **not yet been created or populated**.

See [ORGANIZATION.md](ORGANIZATION.md) for the destination mapping, source revisions, and migration checklist.

For future work, use a permanent project home with a useful default branch, and short-lived task branches within that project.
