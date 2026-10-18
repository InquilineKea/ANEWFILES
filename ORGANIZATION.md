# Project organization and migration

Status: index published; repository split pending. The connected GitHub tools used for this pass support file edits but do not expose repository creation. No destination repository is implied to exist.

## Destination mapping

| Proposed repository | Source branch | Paths to carry forward |
|---|---|---|
| interactive-science-lab | claude/gramian-interactive-simulation-s90ud9 | gramian-lab/ |
| interactive-science-lab | claude/zealous-ptolemy-5ta9rh | noggin-island/ |
| interactive-science-lab | claude/gifted-meitner-cppi8d | computational-matter-zoo/ |
| kagome-research | claude/nice-cerf-umzz6k | analysis/, agsci_tasks/, both twisted-kagome HTML files, applicable .gitignore rules |
| brain-atlas-analysis | claude/endothelial-atlas-download-7vk2je | Atlas scripts, README files, results, and documented dataset provenance |
| boston-environment-analysis | claude/explore-boston-data-mcp-UuLpd | Boston scripts, maps, results, and documented data sources |
| personal-research-tools | claude/check-gtex-eqtl-variants-em0bA | check_gtex_eqtl.py and gtex_eqtl_results.json under gtex-eqtl/ |
| ANEWFILES archive | claude/elegant-noether-tt38nx | Keep existing evaluation branch linked from the index |
| ANEWFILES archive | claude/curl-forum-thread-umkeC | Forum text and project index |

Computational matter zoo is currently a report plus Python checks, so it needs its own Python environment even when grouped with browser demonstrations.

## Source checkpoints

These are the observed branch tips on 2026-10-18, before the index changes. Recheck branch tips before migrating to capture newer work.

| Source | Commit |
|---|---|
| GTEx | 9950c95750f03ce02096339eb945c1c075cd99e2 |
| Forum archive | 0ded36431ac4975dbc886a5513c2a085e6906705 |
| Elegance evaluation | 8531811ac7797332f63ecc78bb519d828d508fd3 |
| Brain / endothelial atlas | 69b8edacfa652bf0d5e859978abf588f12b621cd |
| Boston environment | d10998c876dbfcf7c7b6d5f87d3f6eacb80a4845 |
| Computational matter zoo | df0859a5ac129f0e2cf7d9455f416a4ecefc84aa |
| Gramian Lab | dc48873ca30fba98d38ec3181100687a66effe76 |
| Kagome | ec5a12acb0b45b432a3e92e839c1fb64d2b810bb |
| Noggin Island | 6965ddd48b8e5c88ecccc878c1c2e0a3145fb443 |

## Migration procedure

1. Create the destination repositories with deliberate visibility and a permanent default branch.
2. Refresh and preserve source branch references. Keep ANEWFILES and its historical branches available.
3. Copy the listed project paths with a provenance note linking to the exact source commits. A snapshot copy preserves provenance through the archive but does not transfer full per-file Git history; use history filtering if full history in each destination is required.
4. Keep existing internal paths intact initially. Move scripts, data, and outputs only after checking relative-path assumptions.
5. Audit imports and external tools. Add project-specific environment definitions; produce lockfiles from a successfully resolved and tested environment. Do not invent version pins.
6. Document data sources, download instructions, synthetic fixtures, and generated results. Assess large datasets before copying them into new Git histories.
7. Run each project's appropriate check in its own environment. Record commands, versions, and results in its README; keep checks independent of live API access where feasible.
8. Update the ANEWFILES index to the verified destination URLs. Mark the old branches as historical in documentation without deleting them.

## Project conventions

Each active project should have:

- A README with purpose, example output, setup, run command, and scientific limitations.
- A clearly identified entry point.
- Dependency declarations and a tested lockfile where dependencies are substantial.
- Small reproducible example data, with synthetic examples labeled explicitly.
- An appropriate automated check after its runtime has been verified.
- Separate locations for source, tests, and generated results when moving files will not break existing scripts.
- Source commit and data provenance for imported work.

Branches should describe temporary changes, such as `fix/eigenvalue-validation`. Independent projects should have permanent folders or repositories.

## Validation of this organization pass

The index was built from the nine observed branches and their file trees. Relevant setup files, scripts, and reports were inspected. This pass changes documentation only; no scientific analyses or project runtimes were executed.
