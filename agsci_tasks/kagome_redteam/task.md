# Goal
- Find parameter sets where the Weyl-node detector of the Twisted Kagome–Weyl Hyperlattice page gives the wrong node count.
- Score: one point per distinct verified failure (near-duplicates count once, compared in physically meaningful coordinates: tz is ignored because it never changes a gap, flux is folded with 360−flux, and mA, dz are measured in units of t), plus 0.5 if you find both kinds: missed nodes and false nodes. Higher is better.
- Then explain the failure modes. A good result says *why* the detector fails (which assumption breaks) and proposes a fix, rather than listing random hits.

# Model
- Stacked kagome block, 3×3 Bloch Hamiltonian: nearest-neighbour hopping `-t` with breathing (t▲ = t(1-δ), t▽ = t(1+δ)), phase `flux/3` per bond (flux per triangle in degrees), sublattice mass `mA` on A, and vertical hopping `tz+dz` on A and `tz-dz/2` on B and C, entering as `2·(…)·cos kz`.
- Bounds: t [0.2,2], breath [-0.45,0.45], flux [0,360], mA [-3,3], tz [-1,1], dz [-1.5,1.5]; gap '12' (bands 1|2) or '23' (bands 2|3).

# The detector under test
- The page's current detector: for 49 planes kz in [-pi, pi] it computes the Chern number of the band below the gap (18x18 Fukui grid). It searches for gap minima in windows around Chern jumps (widened by one plane) and in runs of planes whose gap is below 0.01 t. A touching counts only if its refined gap converges below 1e-3 t, and it is a Weyl node only if the Berry flux through a small sphere around it is nonzero.
- History: an earlier version counted Chern jumps directly. An AI-Scientist-v3 red team found its two failure modes (missed nodes behind the gap guard; false nodes from Chern-grid flicker near time reversal). This version fixes both and agrees with the referee on 44/44 validation cases (analysis/weyl_detector_validation.json).

# Referee
- A 56^3 brute-force scan finds band touchings (gap < 1e-3 after refinement); each is classified by its Berry charge through spheres of radius 0.03 and 0.01. Cases where the two radii disagree are "inconclusive" and cannot score. A fine Chern-number count is not used: it was shown to flicker in exactly the hard cases.

# Contract
- `program_0/src/program.py` defines `solution()` returning a list of at most 12 dicts. Runtime limit 900 s.
- `kagome_core.py` beside the program contains `weyl_detect_like_js` and `weyl_gap_scan`; use them to search.
