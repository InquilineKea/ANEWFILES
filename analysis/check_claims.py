"""Re-test the quantitative claims made about the two kagome pages, with the numpy core.

Usage:  python analysis/check_claims.py [--quick]
Writes analysis/claims_report.json. Each check records the claim as stated, what was measured,
and a verdict: holds / partly holds / fails / untested.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from kagome_core import (MoireParams, WeylParams, isolated_band_scores, ldos_window, localisation_metrics,  # noqa: E402
                         manifold_stats, moire_bands, moire_build, moire_grid, moire_path, registry,
                         weyl_detect_like_js, weyl_gap_scan)

QUICK = "--quick" in sys.argv
report = []


def record(cid, claim, measured, verdict, note=""):
    report.append({"id": cid, "claim": claim, "measured": measured, "verdict": verdict, "note": note})
    print(f"[{verdict:>12}] {cid}: {note or claim}")


def dense_ks(sc, nseg=14, nk=8):
    return np.concatenate([moire_path(sc, nseg), moire_grid(sc, nk)])


def narrowest_state_metrics(sc, E):
    """LDOS of the narrowest manifold band at Gamma, window = that band's own range."""
    st = manifold_stats(sc, E)
    b = st["narrowest_band"]
    lo, hi = E[:, b].min(), E[:, b].max()
    E0, win = (lo + hi) / 2, max((hi - lo) / 2, 1e-6) + 1e-9
    wt, ns = ldos_window(sc, sc.high_symmetry()["G"], E0, win)
    return st, {**localisation_metrics(sc, wt), "states": ns, "E0": float(E0)}


t0 = time.time()
# ----------------------------------------------------------------------------- moire: default numbers
sc = moire_build(MoireParams(m=3))
E_page = moire_bands(sc, np.concatenate([moire_path(sc, 14), moire_grid(sc, 4)]))   # what the page samples
E_dense = moire_bands(sc, dense_ks(sc, 20, 10))
st_page, st_dense = manifold_stats(sc, E_page), manifold_stats(sc, E_dense)
record("M1", "m=3 defaults: manifold width 0.425 t, 17 gapped pieces, narrowest band 4.3e-4 t",
       {"page_sampling": st_page, "dense_sampling": st_dense},
       ("holds" if st_dense["pieces"] == 17 else "partly holds") if abs(st_page["manifold_width"] - .425) < .002 and st_page["pieces"] == 17 else "fails",
       "piece count depends on k sampling (sparse sampling overcounts gaps). " f"page sampling: W={st_page['manifold_width']:.4f}, pieces={st_page['pieces']}, narrowest={st_page['narrowest_width']:.2e}; "
       f"dense (20/seg + 10x10): W={st_dense['manifold_width']:.4f}, pieces={st_dense['pieces']}, narrowest={st_dense['narrowest_width']:.2e}")

# ----------------------------------------------------------------------------- AA localisation
_, loc = narrowest_state_metrics(sc, E_dense)
record("M2", "Narrowest manifold states sit on AA patches (x3.4 over an even spread), ~17% participation",
       loc, "holds" if loc["aa_enhancement"] > 2 else "fails",
       f"AA enhancement x{loc['aa_enhancement']:.2f}, participation {100*loc['participation']:.1f}% ({loc['states']} states)")

# ----------------------------------------------------------------------------- twist sweep with dense sampling
sweep = []
for m in range(1, 4 if QUICK else 6):
    s = moire_build(MoireParams(m=m))
    nk = 8 if m <= 3 else 5
    E = moire_bands(s, dense_ks(s, 10, nk))
    st, lc = narrowest_state_metrics(s, E)
    sweep.append({"m": m, "theta_deg": math.degrees(s.theta), "L_M": float(np.hypot(*s.S1)), **st, **{f"loc_{k}": v for k, v in lc.items()}})
page_sweep = {1: .0096, 2: .0019, 3: 2.7e-4, 4: 4.9e-5, 5: 6.3e-6}
ratios = [r["narrowest_width"] / page_sweep[r["m"]] for r in sweep]
LM = np.array([r["L_M"] for r in sweep]); lw = np.log([r["narrowest_width"] for r in sweep])
slope, icpt = np.polyfit(LM, lw, 1)
resid = lw - (slope * LM + icpt)
record("M3", "Narrowest sub-band width falls ~exponentially with moire period (1e-2 t at 21.8 deg to 6e-6 t at 6 deg); 5-k values are lower bounds",
       {"sweep": sweep, "dense_over_5k_ratio": ratios, "fit_log_width_vs_LM": {"slope": slope, "rms_resid": float(np.sqrt(np.mean(resid**2)))}},
       "holds" if slope < 0 and all(r >= 0.99 for r in ratios) else "partly holds",
       "dense/5-k width ratios " + ", ".join(f"m={r['m']}:{x:.2f}" for r, x in zip(sweep, ratios)) +
       f"; log(width) vs L_M slope {slope:.2f}/a, rms residual {np.sqrt(np.mean(resid**2)):.2f}")

# ----------------------------------------------------------------------------- flux delocalises
s60 = moire_build(MoireParams(m=3, flux=60, V=0.2))
mono60 = moire_build(MoireParams(m=1, flux=60, tp=0))
Em60 = moire_bands(mono60, dense_ks(mono60, 8, 6))
mono_top60 = float(Em60[:, -1].max() - Em60[:, -1].min())
E60 = moire_bands(s60, dense_ks(s60, 14, 6))
st60, lc60 = narrowest_state_metrics(s60, E60)
record("M4", "Flux 60 deg (V=0.2) delocalises the flat states and kills the AA preference (manifold ~1.45 t)",
       {"stats": st60, "loc": lc60, "monolayer_top_band_width_at_flux60": mono_top60},
       "holds" if lc60["aa_enhancement"] < 1.5 else "fails",
       f"narrowest band now: AA x{lc60['aa_enhancement']:.2f}, participation {100*lc60['participation']:.1f}%, manifold {st60['manifold_width']:.3f} t. "
       "The page's x0.92 came from a stale LDOS energy (bug since fixed). Flux does widen the narrowest band ~10x. "
       f"The folded monolayer top band is no longer flat ({mono_top60:.2f} t wide, measured in the m=1 cell), yet AA confinement survives.")

# ----------------------------------------------------------------------------- the discriminating experiment
rows = []
for delta in ([0, .2] if QUICK else [0, .1, .2, .3]):
    for tp in ([.35] if QUICK else [.15, .35, .6]):
        s = moire_build(MoireParams(m=3, breath=delta, tp=tp))
        E = moire_bands(s, dense_ks(s, 10, 6))
        st, lc = narrowest_state_metrics(s, E)
        mono = moire_build(MoireParams(m=1, breath=delta, tp=0))
        Em = moire_bands(mono, dense_ks(mono, 8, 6))
        top = Em[:, -1]
        rows.append({"delta": delta, "tp": tp, "aa": lc["aa_enhancement"], "part": lc["participation"],
                     "narrowest": st["narrowest_width"], "monolayer_flat_width": float(top.max() - top.min())})
aa_persist = all(r["aa"] > 1.5 for r in rows)
record("M5", "Untested before: does the AA preference survive breathing (delta != 0, flux 0)? CLS reading predicts yes",
       rows, "holds" if aa_persist else "partly holds",
       "; ".join(f"d={r['delta']},t⊥={r['tp']}: AA x{r['aa']:.2f}" for r in rows) +
       f". Monolayer flat band under breathing stays flat (width {max(r['monolayer_flat_width'] for r in rows):.1e} t), so delta does not distinguish the two readings as cleanly as I implied.")

# ----------------------------------------------------------------------------- Weyl detector vs brute force
wp = WeylParams()
det12 = weyl_detect_like_js(wp, "12")
det23 = weyl_detect_like_js(wp, "23")
bf12 = [n for n in weyl_gap_scan(wp, 0, 40) if n[3] < 1e-4]
bf23 = [n for n in weyl_gap_scan(wp, 1, 40) if n[3] < 1e-4]
record("W1", "Defaults (flux 90, dz 0.8): 2 nodes on gap 1|2 at kz = ±2.556, 2 nodes on gap 2|3 at kz = ±0.59",
       {"detector_12": det12["n_nodes"], "detector_23": det23["n_nodes"], "brute_12": bf12, "brute_23": bf23},
       "holds" if det12["n_nodes"] == len(bf12) == 2 and det23["n_nodes"] == len(bf23) == 2 else "fails",
       f"detector {det12['n_nodes']}/{det23['n_nodes']}, brute force {len(bf12)}/{len(bf23)}; brute kz: "
       + ", ".join(f"{n[2]:+.3f}" for n in bf12 + bf23))

rng = np.random.default_rng(1)
trials = []
for _ in range(10 if QUICK else 40):
    p = WeylParams(t=1, breath=float(rng.uniform(-.3, .3)), flux=float(rng.uniform(0, 360)), mA=float(rng.uniform(-2, 2)),
                   tz=float(rng.uniform(-.6, .6)), dz=float(rng.uniform(-1.5, 1.5)))
    d = weyl_detect_like_js(p, "12")
    b = [n for n in weyl_gap_scan(p, 0, 32) if n[3] < 1e-3]
    trials.append({"params": p.__dict__, "detector": d["n_nodes"], "brute": len(b), "ungapped_planes": d["ungapped_planes"]})
# referee every disagreement at higher resolution (finer Chern grid and kz planes, finer 3D scan)
for t in trials:
    if t["detector"] != t["brute"]:
        p = WeylParams(**t["params"])
        t["referee_brute"] = len([n for n in weyl_gap_scan(p, 0, 56) if n[3] < 1e-3])
        t["referee_detector_fine"] = weyl_detect_like_js(p, "12", nk=30, nz=96)["n_nodes"]
        t["fault"] = "brute-force scan (too coarse)" if t["referee_brute"] == t["detector"] else "detector"
agree = sum(t["detector"] == t["brute"] for t in trials)
det_right = sum(t["detector"] == t.get("referee_brute", t["brute"]) for t in trials)
miss = [t for t in trials if t["detector"] < t["brute"]]
extra = [t for t in trials if t["detector"] > t["brute"]]
record("W2", "Implied: the C(kz)-jump detector finds the Weyl nodes (no stated error rate)",
       {"agree_first_pass": agree, "detector_correct_after_referee": det_right, "n": len(trials), "missed": miss, "extra": extra},
       "holds" if det_right == len(trials) else "partly holds",
       f"first pass {agree}/{len(trials)} agree; after refereeing disagreements at higher resolution the page's detector is right in "
       f"{det_right}/{len(trials)}. Faults: " + ", ".join(t["fault"] for t in miss + extra))

# ----------------------------------------------------------------------------- analytic claims
record("S1", "Twist-locked interlayer maps (kappa) are pure gauge: H0 unchanged",
       "omega = kappa*(alpha_l - alpha_{l+1}) = g_{l+1} - g_l with g_l = -kappa*alpha_l, so every cycle holonomy cancels",
       "holds", "exact by construction (a coboundary); no numerics needed")
record("S2", "Spin gauge with gamma = 1 is a pure Z2 gauge",
       "Rot(omega + pi(1 - s_i s_j)/2) = s_i s_j Rot(omega); x_i -> s_i x_i removes it",
       "holds", "exact; fractional gamma gives holonomy pi*gamma*(#anti-aligned edges), which can be nonzero")

out = Path(__file__).parent / "claims_report.json"
out.write_text(json.dumps({"elapsed_s": time.time() - t0, "quick": QUICK, "checks": report}, indent=1, default=float))
print(f"\nwrote {out} in {time.time()-t0:.0f} s")
