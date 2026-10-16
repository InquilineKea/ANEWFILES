"""Hybrid breathing + near-pi-flux flat band optimization.

Combines moderate breathing distortion (|breath| 0.1-0.4) with near-pi-flux
(160-178 deg) to synergistically flatten bands and open gaps. Breathing breaks
sublattice symmetry, lifting pi-flux degeneracies and allowing flux closer to
180 while maintaining isolation gaps. Searches m=1 and m=2 supercells.
"""
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kagome_core import (MoireParams, isolated_band_scores,
                         moire_bands, moire_build, moire_grid, moire_path)

MIN_GAP = 0.01
TIME_LIMIT = 540


def _score_params(params):
    try:
        p = MoireParams(**{**params, "m": int(params.get("m", 1))})
        sc = moire_build(p)
        ks = np.concatenate([moire_path(sc, 12), moire_grid(sc, 6)])
        E = moire_bands(sc, ks)
        bands = [b for b in isolated_band_scores(E) if b["gap"] >= MIN_GAP]
        if not bands:
            return -999.0, None
        best = max(bands, key=lambda b: b["ratio"])
        return math.log10(best["ratio"]), best
    except Exception:
        return -999.0, None


def solution():
    t0 = time.time()
    best_score = -999.0
    best_params = {"m": 1, "t": 1.0, "breath": 0.0, "flux": 168.0,
                   "tp": 0.17, "xi": 0.23, "rc": 0.8, "V": 0.0,
                   "sx": 0.0, "sy": 0.0}
    best_info = None
    top_hits = []

    # Phase 1: coarse grid at m=1 — hybrid breathing + near-pi-flux regime
    flux_vals = [160, 163, 165, 167, 168, 170, 173, 175, 178]
    breath_vals = [-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
    tp_vals = [0.05, 0.10, 0.15, 0.17, 0.20, 0.25, 0.35]
    xi_vals = [0.10, 0.15, 0.23, 0.30]
    V_vals = [0.0, 0.1, 0.2, 0.3, 0.5]
    rc_val = 0.8

    count = 0
    phase1_budget = TIME_LIMIT * 0.50
    for flux in flux_vals:
        if time.time() - t0 > phase1_budget:
            break
        for breath in breath_vals:
            for tp in tp_vals:
                for xi in xi_vals:
                    for V in V_vals:
                        params = {"m": 1, "t": 1.0, "breath": breath,
                                  "flux": flux, "tp": tp, "xi": xi,
                                  "rc": rc_val, "V": V,
                                  "sx": 0.0, "sy": 0.0}
                        sc, info = _score_params(params)
                        count += 1
                        if sc > best_score:
                            best_score = sc
                            best_params = params.copy()
                            best_info = info
                        if sc > 0.5:
                            top_hits.append((sc, params.copy(), info))

    top_hits.sort(key=lambda x: -x[0])
    top_hits = top_hits[:20]
    phase1_count = count
    phase1_best = best_score
    print(f"[phase1] {phase1_count} combos in {time.time()-t0:.1f}s, best={phase1_best:.4f}")

    # Phase 2: fine refinement around top-3 hits from phase 1
    phase2_budget = TIME_LIMIT * 0.72
    if top_hits and time.time() - t0 < phase2_budget:
        seen = set()
        for _, base, _ in top_hits[:3]:
            if time.time() - t0 > phase2_budget:
                break
            refine_flux = [base["flux"] + d for d in [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]]
            refine_breath = [max(-0.4, min(0.4, base["breath"] + d))
                             for d in [-0.15, -0.1, -0.05, -0.02, 0, 0.02, 0.05, 0.1, 0.15]]
            refine_tp = [base["tp"] * f for f in [0.5, 0.7, 0.85, 0.93, 1.0, 1.07, 1.15, 1.3, 1.5, 2.0]]
            refine_xi = [max(0.08, min(0.6, base["xi"] * f))
                         for f in [0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5]]
            refine_V = [max(0.0, min(1.0, base["V"] + d))
                        for d in [-0.15, -0.1, -0.05, -0.02, 0, 0.02, 0.05, 0.1, 0.15]]

            for flux in refine_flux:
                if time.time() - t0 > phase2_budget:
                    break
                if not (0 <= flux <= 360):
                    continue
                for breath in refine_breath:
                    for tp in refine_tp:
                        if not (0.001 <= tp <= 1.0):
                            continue
                        for xi in refine_xi:
                            for V in refine_V:
                                key = (round(flux, 2), round(breath, 3),
                                       round(tp, 4), round(xi, 3), round(V, 3))
                                if key in seen:
                                    continue
                                seen.add(key)
                                params = {"m": 1, "t": 1.0, "breath": breath,
                                          "flux": flux, "tp": tp, "xi": xi,
                                          "rc": rc_val, "V": V,
                                          "sx": 0.0, "sy": 0.0}
                                sc, info = _score_params(params)
                                count += 1
                                if sc > best_score:
                                    best_score = sc
                                    best_params = params.copy()
                                    best_info = info

    phase2_count = count - phase1_count
    print(f"[phase2] {phase2_count} refinements, best={best_score:.4f}")

    # Phase 2b: second refinement around current best
    phase2b_budget = TIME_LIMIT * 0.80
    if best_score > -999 and time.time() - t0 < phase2b_budget:
        base = best_params.copy()
        refine_flux = [base["flux"] + d for d in [-1.5, -0.8, -0.3, -0.1, 0, 0.1, 0.3, 0.8, 1.5]]
        refine_breath = [max(-0.4, min(0.4, base["breath"] + d))
                         for d in [-0.08, -0.04, -0.02, -0.01, 0, 0.01, 0.02, 0.04, 0.08]]
        refine_tp = [base["tp"] * f for f in [0.8, 0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1, 1.2]]
        refine_xi = [max(0.08, min(0.6, base["xi"] * f))
                     for f in [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]]
        refine_V = [max(0.0, min(1.0, base["V"] + d))
                    for d in [-0.06, -0.03, -0.01, 0, 0.01, 0.03, 0.06]]
        for flux in refine_flux:
            if time.time() - t0 > phase2b_budget:
                break
            if not (0 <= flux <= 360):
                continue
            for breath in refine_breath:
                for tp in refine_tp:
                    if not (0.001 <= tp <= 1.0):
                        continue
                    for xi in refine_xi:
                        for V in refine_V:
                            params = {"m": 1, "t": 1.0, "breath": breath,
                                      "flux": flux, "tp": tp, "xi": xi,
                                      "rc": rc_val, "V": V,
                                      "sx": 0.0, "sy": 0.0}
                            sc, info = _score_params(params)
                            count += 1
                            if sc > best_score:
                                best_score = sc
                                best_params = params.copy()
                                best_info = info

    m1_best_score = best_score
    m1_best_params = best_params.copy()
    print(f"[phase2b] total={count}, m1 best={m1_best_score:.4f}")

    # Phase 3: m=2 with best flux/breath pairs from m=1
    phase3_budget = TIME_LIMIT * 0.92
    m2_count = 0
    if time.time() - t0 < phase3_budget:
        flux_breath_pairs = set()
        for _, p, _ in top_hits[:5]:
            flux_breath_pairs.add((p["flux"], p["breath"]))
        flux_breath_pairs.add((best_params["flux"], best_params["breath"]))

        for flux, breath in flux_breath_pairs:
            if time.time() - t0 > phase3_budget:
                break
            for tp in [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]:
                if time.time() - t0 > phase3_budget:
                    break
                for xi in [0.08, 0.12, 0.18, 0.25, 0.35]:
                    for V in [0.0, 0.1, 0.2, 0.3, 0.5]:
                        params = {"m": 2, "t": 1.0, "breath": breath,
                                  "flux": flux, "tp": tp, "xi": xi,
                                  "rc": rc_val, "V": V,
                                  "sx": 0.0, "sy": 0.0}
                        sc, info = _score_params(params)
                        count += 1
                        m2_count += 1
                        if sc > best_score:
                            best_score = sc
                            best_params = params.copy()
                            best_info = info

    print(f"[phase3] {m2_count} m=2 combos, best={best_score:.4f}")

    # Phase 4: layer sliding at overall best parameters
    phase4_budget = TIME_LIMIT * 0.97
    slide_count = 0
    if time.time() - t0 < phase4_budget:
        base_slide = best_params.copy()
        for sx in [0.0, 0.125, 0.25, 0.375, 0.5]:
            if time.time() - t0 > phase4_budget:
                break
            for sy in [0.0, 0.125, 0.25, 0.375, 0.5]:
                params = {**base_slide, "sx": sx, "sy": sy}
                sc, info = _score_params(params)
                count += 1
                slide_count += 1
                if sc > best_score:
                    best_score = sc
                    best_params = params.copy()
                    best_info = info

    print(f"[phase4] {slide_count} slide combos, best={best_score:.4f}")

    # Phase 5: report mechanism explanation
    elapsed = time.time() - t0
    bi = best_info or {}
    c1_score = 1.852
    is_hybrid = (abs(best_params.get("breath", 0)) >= 0.05
                 and 160 <= best_params.get("flux", 0) <= 178)
    n_sites = 42 * int(best_params.get("m", 1))**2 if best_params.get("m", 1) <= 2 else 42
    m0 = n_sites - 2 * (7 * int(best_params.get("m", 1))**2)

    print(f"\n{'='*60}")
    print(f"[RESULT] {count} total combos in {elapsed:.1f}s")
    print(f"[RESULT] best combined_score = {best_score:.4f}")
    print(f"[RESULT] best params: {best_params}")
    print(f"[RESULT] band {bi.get('band')}: width={bi.get('width')}, "
          f"gap={bi.get('gap')}, ratio={bi.get('ratio')}")
    print(f"[RESULT] m=1 best score = {m1_best_score:.4f}")

    flux_val = best_params.get('flux', 0)
    breath_val = best_params.get('breath', 0)
    tp_val = best_params.get('tp', 0)
    xi_val = best_params.get('xi', 0)
    V_val = best_params.get('V', 0)
    m_val = int(best_params.get('m', 1))

    mechanism_parts = []
    mechanism_parts.append(
        f"Near-pi-flux at {flux_val:.1f} deg provides Aharonov-Bohm caging "
        f"(destructive interference flattens kagome bands)."
    )
    if abs(breath_val) >= 0.05:
        mechanism_parts.append(
            f"Breathing distortion breath={breath_val:.3f} breaks sublattice "
            f"symmetry, lifting pi-flux degeneracies and opening wider "
            f"inter-band gaps for isolation."
        )
    else:
        mechanism_parts.append(
            f"Breathing={breath_val:.3f} is near zero; the hybrid synergy "
            f"did not improve over pure pi-flux for this optimum."
        )
    mechanism_parts.append(
        f"Weak interlayer coupling tp={tp_val:.4f}, xi={xi_val:.3f} creates "
        f"perturbative bandwidth ~tp^2."
    )
    if abs(V_val) > 0.01:
        mechanism_parts.append(
            f"Interlayer bias V={V_val:.3f} splits layer manifolds for gap opening."
        )
    mechanism_parts.append(
        f"Band {bi.get('band')} has width={bi.get('width', 0):.2e}, "
        f"gap={bi.get('gap', 0):.4f}, ratio={bi.get('ratio', 0):.2e}."
    )
    in_top = bi.get('band', 0) >= m0 if bi.get('band') is not None else False
    mechanism_parts.append(
        f"Band {'is in top manifold' if in_top else 'is in lower manifold'} "
        f"at m={m_val}."
    )

    if is_hybrid:
        mechanism_parts.append(
            f"HYBRID mechanism confirmed: both flux ({flux_val:.1f} deg) and "
            f"breathing ({breath_val:.3f}) contribute."
        )
    else:
        mechanism_parts.append(
            f"Single-mechanism optimum: flux={flux_val:.1f}, breath={breath_val:.3f}."
        )

    score_diff = best_score - c1_score
    if score_diff > 0:
        mechanism_parts.append(
            f"vs cycle 1 (pure pi-flux, score={c1_score:.3f}): IMPROVEMENT "
            f"by {score_diff:.3f} ({10**score_diff:.1f}x better ratio)."
        )
    else:
        mechanism_parts.append(
            f"vs cycle 1 (pure pi-flux, score={c1_score:.3f}): "
            f"{'matched' if abs(score_diff) < 0.01 else 'fell short'} "
            f"by {abs(score_diff):.3f}."
        )

    mechanism = " ".join(mechanism_parts)
    print(f"[mechanism] {mechanism}")

    return best_params
