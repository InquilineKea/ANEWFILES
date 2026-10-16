# Goal
- Find parameters of a commensurate twisted kagome bilayer that produce the flattest *isolated* electronic band.
- Score: `combined_score = log10(gap / width)` for the best band whose energy gap to **both** neighbouring bands is at least 0.01 t. Higher is better.
- Explore genuinely different mechanisms, not one knob tuned harder. Candidate routes include moiré confinement in AA patches at small twist, the kagome compact-localised-state flat band, flux-induced Chern bands, breathing distortion, interlayer bias, layer sliding, and the range of the interlayer hopping.

# Model
- Two kagome sheets with nearest-neighbour hopping `-t` (t = 1 fixed), rotated by ∓θ/2, at commensurate angles `cos θ = (3m²+3m+½)/(3m²+3m+1)`, m = 1..4 (21.8°, 13.2°, 9.43°, 7.34°).
- Interlayer hopping `tp · exp(-ρ²/2ξ²)` for in-plane distance ρ < min(rc, 4ξ).
- Parameters and bounds: `m` int [1,4]; `breath` [-0.4,0.4] (t▲ = t(1-δ), t▽ = t(1+δ)); `flux` [0,360] degrees per triangle; `tp` [-1,1]; `xi` [0.08,0.6]; `rc` [0.2,1.2]; `V` [-1,1] interlayer bias; `sx`, `sy` [0,1] slide of layer 1.

# Contract
- `program_0/src/program.py` defines `solution()` returning a dict with any subset of the parameters above (missing ones take their defaults).
- The evaluator re-solves the Bloch Hamiltonian on a Γ–K–M–Γ path (12 points per segment) plus a 6×6 moiré-zone grid. Widths and gaps are therefore measured on that sampling; bands are compared by index after sorting.
- `kagome_core.py` next to the program is a copy of the solver you may use or modify for searching. Runtime limit 600 s.

# What to report
- Alongside the score, explain *why* the winning band is flat: does it live on AA patches (`aa_enhancement` > 1), is it a descendant of the monolayer flat band (`best_band_in_top_manifold`), and does it survive switching flux on or off? A good result is one whose mechanism is understood, not only one with a high ratio.
