# Claim checks

`kagome_core.py` is an independent numpy implementation of both pages' models. It agrees with
the moiré page's JavaScript solver to 2e-14 on eigenvalues across bias, flux, breathing, slide
and negative t⊥. `check_claims.py` re-tests the quantitative claims made about the pages and
writes `claims_report.json` (about 30 s; add `--quick` for a smaller run). Needs numpy only.

Results from the last run:

| ID | Claim | Verdict |
|---|---|---|
| M1 | m=3 defaults: manifold 0.425 t, 17 gapped pieces, narrowest 4.3e-4 t | partly holds: width reproduces; denser sampling gives 13 pieces and 4.7e-4 t. Piece counts depend on k sampling |
| M2 | Narrowest manifold states sit on AA patches | holds: ×6.8 over an even spread (the page's ×3.4 used a wider window holding 2 states) |
| M3 | Narrowest width falls ~exponentially with moiré period | holds: slope −1.11 per a of L_M. The page's 5-k values are lower bounds, 1.1–1.8× too small |
| M4 | Flux 60° delocalises the flat states and removes the AA preference | **fails**: AA ×4.4 survives. The earlier ×0.92 came from a stale-energy bug, since fixed |
| M5 | (new) AA preference survives breathing δ | holds (×5.5–6.8), but breathing keeps the monolayer flat band exactly flat, so it does not discriminate |
| W1 | Default Weyl block: 2 + 2 nodes at kz = ±2.556, ±0.586 | holds (brute-force 3D scan agrees) |
| W2 | The Chern-jump detector finds the nodes | partly holds: right in 37/40 random cases; 2 misses (gap guard), 1 false positive (coarse grid) |
| S1, S2 | κ twist-lock and γ = 1 spin gauge are pure gauge | hold, exactly (analytic) |

M4 plus M5 shift the interpretation. AA confinement survives flux, which destroys the monolayer
flat band, so the narrow states are not just compact localised states selected by stacking.
Moiré-confined bound states in the AA potential fit better, together with the exponential
width-vs-L_M trend.
