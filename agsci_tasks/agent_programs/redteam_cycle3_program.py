"""False-node priority search with slot reservation for Weyl detector failures.

Searches for false-node overcounts FIRST (page > referee), reserving up to 3
of 12 output slots exclusively for them, then fills remaining slots with
missed-node failures from the proven tz-sweep baseline.

False-node mechanism: At small magnetic flux (1-15 deg), Berry curvature
concentrates in a BZ region smaller than the 18x18 FHS grid spacing
(~0.35 rad). The coarse discretization yields wrong Chern integers, creating
spurious transitions the 30x30 referee resolves correctly. Large |mA| and
nonzero breath further distort the band structure, amplifying this effect.
Fix: use >= 24x24 Chern grid or adaptive mesh refinement.

Missed-node mechanism: A Weyl node near a sampled kz plane makes the band gap
< gap_tol=0.01, causing the detector to filter out flanking Chern transitions.
Fix: adaptive gap threshold or interpolate Chern across gapped planes.
"""

import time

import numpy as np

from kagome_core import WeylParams, weyl_detect_like_js

SCALE = np.array([1.8, 0.9, 360.0, 6.0, 2.0, 3.0])
DEDUP = 0.03
TIME_BUDGET = 600
MAX_RESULTS = 12
FALSE_SLOT_CAP = 3


def _vec(t, breath, flux, mA, tz, dz):
    return np.array([t, breath, flux, mA, tz, dz]) / SCALE


def _distinct(v, gap, found):
    return all(np.linalg.norm(v - f["v"]) > DEDUP for f in found if f["gap"] == gap)


def _check(t, breath, flux, mA, tz, dz, gap, found, t0):
    if time.time() - t0 > TIME_BUDGET:
        return None
    v = _vec(t, breath, flux, mA, tz, dz)
    if not _distinct(v, gap, found):
        return None
    try:
        p = WeylParams(t=t, breath=breath, flux=flux, mA=mA, tz=tz, dz=dz)
        page = weyl_detect_like_js(p, gap)["n_nodes"]
        if time.time() - t0 > TIME_BUDGET:
            return None
        ref = weyl_detect_like_js(p, gap, nk=30, nz=96)["n_nodes"]
        if page != ref:
            return {
                "params": {"t": t, "breath": breath, "flux": flux,
                           "mA": mA, "tz": tz, "dz": dz, "gap": gap},
                "v": v, "gap": gap,
                "kind": "false" if page > ref else "miss",
            }
    except Exception:
        pass
    return None


def solution():
    t0 = time.time()
    false_found = []
    miss_found = []

    def elapsed():
        return time.time() - t0

    def all_found():
        return false_found + miss_found

    def try_false(t, breath, flux, mA, tz, dz, gap):
        if elapsed() > TIME_BUDGET or len(false_found) >= FALSE_SLOT_CAP:
            return
        r = _check(t, breath, flux, mA, tz, dz, gap, all_found(), t0)
        if r and r["kind"] == "false":
            false_found.append(r)
        elif r and r["kind"] == "miss":
            miss_cap = MAX_RESULTS - max(len(false_found), 1)
            if len(miss_found) < miss_cap:
                miss_found.append(r)

    def try_miss(t, breath, flux, mA, tz, dz, gap):
        miss_cap = MAX_RESULTS - len(false_found)
        if elapsed() > TIME_BUDGET or len(miss_found) >= miss_cap:
            return
        r = _check(t, breath, flux, mA, tz, dz, gap, all_found(), t0)
        if r and r["kind"] == "false":
            if len(false_found) < FALSE_SLOT_CAP:
                false_found.append(r)
        elif r and r["kind"] == "miss":
            if len(miss_found) < miss_cap:
                miss_found.append(r)

    # === PHASE 1: FALSE-NODE SEARCH (runs first, budget ~400s) ===
    false_budget = 400

    # Phase 1a: Small flux (1-15 deg) with large |mA|
    small_fluxes = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15]
    large_mAs = [-3.0, -2.5, -2.0, -1.5, 1.5, 2.0, 2.5, 3.0]
    dzs_false = [0.3, 0.5, 0.8, 1.0]
    breaths_false = [0.0, 0.15, 0.3]
    for flux in small_fluxes:
        for mA in large_mAs:
            for dz in dzs_false:
                for breath in breaths_false:
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                    try_false(1.0, breath, float(flux), mA, 0.0, dz, "12")
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break
        if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
            break

    # Phase 1b: Near-360 flux (time-reversed small flux)
    if len(false_found) < FALSE_SLOT_CAP and elapsed() < false_budget:
        near360_fluxes = [345, 350, 353, 355, 357, 358, 359]
        for flux in near360_fluxes:
            for mA in large_mAs:
                for dz in dzs_false:
                    for breath in breaths_false:
                        if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                            break
                        try_false(1.0, breath, float(flux), mA, 0.0, dz, "12")
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break

    # Phase 1c: Very small dz — kz aliasing with moderate flux
    if len(false_found) < FALSE_SLOT_CAP and elapsed() < false_budget:
        small_dzs = [0.02, 0.05, 0.08, 0.1, 0.15]
        mod_fluxes = [30, 60, 90, 120, 180]
        for dz in small_dzs:
            for flux in mod_fluxes:
                for mA in [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0]:
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                    try_false(1.0, 0.0, float(flux), mA, 0.0, dz, "12")
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break

    # Phase 1d: gap='23' with small flux and large |mA|
    if len(false_found) < FALSE_SLOT_CAP and elapsed() < false_budget:
        for flux in [3, 5, 7, 10, 15]:
            for mA in [-3.0, -2.0, 2.0, 3.0]:
                for dz in [0.5, 0.8, 1.0]:
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                    try_false(1.0, 0.0, float(flux), mA, 0.0, dz, "23")
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break

    # Phase 1e: Small flux with nonzero tz
    if len(false_found) < FALSE_SLOT_CAP and elapsed() < false_budget:
        for flux in [3, 5, 7, 10]:
            for mA in [-3.0, -2.0, 2.0, 3.0]:
                for tz in [-0.5, -0.3, 0.3, 0.5]:
                    for dz in [0.5, 1.0]:
                        if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                            break
                        try_false(1.0, 0.0, float(flux), mA, float(tz), dz, "12")
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break

    # Phase 1f: Moderate flux with extreme params
    if len(false_found) < FALSE_SLOT_CAP and elapsed() < false_budget:
        for flux in [20, 25, 30, 40, 50, 340, 335, 330, 320]:
            for mA in [-3.0, 3.0]:
                for breath in [0.3, 0.4, -0.3, -0.4]:
                    for dz in [0.3, 0.5, 1.0, 1.4]:
                        if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                            break
                        try_false(1.0, breath, float(flux), mA, 0.0, dz, "12")
                    if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                        break
                if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                    break
            if len(false_found) >= FALSE_SLOT_CAP or elapsed() > false_budget:
                break

    # === PHASE 2: MISSED-NODE SEARCH (tz sweep, fills remaining slots) ===
    BT = 0.08479690148362495
    BF = 306.94782185303643
    BM = 0.37176407241713605
    BD = 1.0196445630942264

    # Phase 2a: Known baseline + time-reversed
    try_miss(1.0, BT, BF, BM, -0.28788306271533215, BD, "12")
    try_miss(1.0, BT, 360.0 - BF, BM, -0.28788306271533215, BD, "12")

    # Phase 2b: tz sweep around baseline
    for tz in np.linspace(-0.95, 0.95, 30):
        miss_cap = MAX_RESULTS - len(false_found)
        if elapsed() > TIME_BUDGET or len(miss_found) >= miss_cap:
            break
        try_miss(1.0, BT, BF, BM, float(tz), BD, "12")

    # Phase 2c: Time-reversed tz sweep
    for tz in np.linspace(-0.95, 0.95, 20):
        miss_cap = MAX_RESULTS - len(false_found)
        if elapsed() > TIME_BUDGET or len(miss_found) >= miss_cap:
            break
        try_miss(1.0, BT, 360.0 - BF, BM, float(tz), BD, "12")

    # Phase 2d: gap='23' tz sweep
    for tz in np.linspace(-0.95, 0.95, 15):
        miss_cap = MAX_RESULTS - len(false_found)
        if elapsed() > TIME_BUDGET or len(miss_found) >= miss_cap:
            break
        try_miss(1.0, BT, BF, BM, float(tz), BD, "23")

    # Combine: false nodes first, then missed nodes
    combined = false_found + miss_found
    if not combined:
        return [{"t": 1.0, "breath": BT, "flux": BF,
                 "mA": BM, "tz": -0.28788306271533215, "dz": BD, "gap": "12"}]

    return [f["params"] for f in combined[:MAX_RESULTS]]
