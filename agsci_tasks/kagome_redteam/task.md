# Goal
- Find parameter sets where the Weyl-node detector of the Twisted Kagome–Weyl Hyperlattice page gives the wrong node count.
- Score: one point per distinct verified failure (near-duplicates count once, compared in physically meaningful coordinates: tz is ignored because it never changes a gap, flux is folded with 360−flux, and mA, dz are measured in units of t), plus 0.5 if you find both kinds: missed nodes and false nodes. Higher is better.
- Then explain the failure modes. A good result says *why* the detector fails (which assumption breaks) and proposes a fix, rather than listing random hits.

# Model
- Stacked kagome block, 3×3 Bloch Hamiltonian: nearest-neighbour hopping `-t` with breathing (t▲ = t(1-δ), t▽ = t(1+δ)), phase `flux/3` per bond (flux per triangle in degrees), sublattice mass `mA` on A, and vertical hopping `tz+dz` on A and `tz-dz/2` on B and C, entering as `2·(…)·cos kz`.
- Bounds: t [0.2,2], breath [-0.45,0.45], flux [0,360], mA [-3,3], tz [-1,1], dz [-1.5,1.5]; gap '12' (bands 1|2) or '23' (bands 2|3).

# The detector under test
- For 49 planes kz ∈ [-π, π], compute the Chern number of the band below the gap on an 18×18 grid (Fukui–Hatsugai–Suzuki). Count |ΔC| at each jump between neighbouring planes, but skip jumps where either plane has a minimum gap below 0.01 t.
- Known first-pass behaviour: on 40 random parameter sets it agreed with a careful reference in 37. The 3 failures were 2 misses (gap nearly closed on the neighbouring planes) and 1 false positive (coarse Chern grid).

# Referee
- Each case is checked by (a) the same detector at 30×30 and 97 planes and (b) a 56³ brute-force scan for gap minima below 1e-3 after local refinement. Only cases where (a) and (b) agree can score. If they disagree the case is "inconclusive", which is also worth reporting.

# Contract
- `program_0/src/program.py` defines `solution()` returning a list of at most 12 dicts. Runtime limit 900 s.
- `kagome_core.py` beside the program contains `weyl_detect_like_js` and `weyl_gap_scan`; use them to search.
