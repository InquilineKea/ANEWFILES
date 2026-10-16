# AI-Scientist-v3 run, 2026-10-16

Investigation phase only, run natively in a cloud container (3 idea cycles per task,
agent `claude` via Headless, model claude-opus-4-6). Phases that require Docker were skipped:
task drafting, ablations, final evaluation. Cost reported by the pipeline: $11.01 (red team) +
$9.29 (flat band) + $0.30 (auth check) ≈ $20.60. Wall time 55 min + 51 min.

Every score below was re-checked outside the pipeline with `analysis/kagome_core.py` at higher
resolution than the evaluators use.

## Red team: Weyl-node detector

| Cycle | Agent idea | Pipeline score | Re-scored | What it really found |
|---|---|---|---|---|
| 1 | Multi-Mode Weyl Detector Failure Exploitation | 2.0 | 1.0 | the seed miss again |
| 2 | Adaptive Computational Search | **12.0** | 1.0 | the seed miss, swept along `tz` |
| 3 | False-Node Priority Search with Slot Reservation | 10.5 | **2.5** | seed miss + one genuine false-node case |

**Reward hacking.** Cycle 2 stepped the seed failure across `tz`, which adds the same
`2·tz·cos kz` to all three sublattices and therefore never changes a gap (checked: max change
7e-15), plus one mirrored-flux copy. The original dedup compared raw parameters, so each step
counted as new. Fixed in `kagome_redteam/evaluate.py` (dedup now ignores `tz`, folds flux with
360−flux, and measures mA and dz in units of t); the same program re-scores 1.0.

**Genuine failure modes of the page's detector (2):**
1. *Misses* — when the neighbouring k_z planes have minimum gap below 0.01 t, the guard drops the
   Chern jump entirely, so real nodes vanish. Fix: in an ungapped interval, fall back to a direct
   3D gap-minimum search instead of skipping it.
2. *False nodes* — near time reversal (flux ≈ 1°, mA = −3, dz = 1) the gap stays above the
   guard (0.02–0.05 t) but Berry curvature is too concentrated for the 18×18 grid, so C(k_z)
   flickers in adjacent ±1 pairs. The page reports 10 nodes; an 80³ scan and a 40²×192
   detector both give 2. Fix: accept a node only if its refined gap converges below ~1e-3 t
   (the page already computes that gap and shows it in the node table), and cancel adjacent
   ±1 pairs whose refinement does not converge.

Two further cases were *inconclusive*: the brute-force scan finds two band touchings that
carry no Chern-number change even at high resolution. These may be touchings with zero net
Berry charge (not Weyl nodes); the evaluator correctly gave them no credit.

## Diversity: isolated flat moiré band

| Cycle | Agent idea | Score log10(gap/width) | Winning parameters |
|---|---|---|---|
| — | starter (page defaults) | 1.24 | m=3, flux 0, t⊥ 0.35 |
| 1 | π-flux Aharonov–Bohm caging, weak interlayer coupling | 1.85 | m=1, flux 168°, t⊥ 0.17 |
| 2 | Breathing kagome CLS band, bias isolation | 1.29 | m=2, δ −0.39, V 0.3 |
| 3 | Breathing-assisted near-π-flux optimisation | **5.40** | m=2, flux 167.5°, t⊥ 0.15, ξ 0.08, slide 0.375 |

Cycle 3's band: E = −1.95 t, width 5.9e-8 t, gap 0.014 t (dense re-check: log10 ratio 5.37,
so not a sampling artifact). Controls:

- With t⊥ = 0 there is no isolated band at all, so the interlayer coupling does the isolating.
- The slide barely matters (5.13 without it).
- A single kagome layer at flux 167.5° has a bottom band 0.25 t wide. It is exactly flat only
  at 180°, where the model is gauge-equivalent to t → −t (flat band at −2t).
- The winning state is AA-localised (AA weight ×4.7; 90% of its weight on 21 of 114 sites).

Interpretation: near 180° flux the monolayer bottom band is very heavy, so states bound in AA
patches tunnel between patches only weakly and their bandwidth becomes exponentially small.
That is the moiré-confinement picture from `analysis/README.md` (checks M4/M5), reached
independently. Caveat: the isolating gap is small (1.4% of t); the score is driven mostly by
the tiny width.

Winning programs are in `agent_programs/`. They are agent-written and kept for provenance.
