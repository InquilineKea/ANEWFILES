"""Evaluator: red-team the Weyl-node detector of the Twisted Kagome-Weyl Hyperlattice page.

The detector under test is the page's current one (weyl_detect_v2): Chern-jump windows plus
nearly-gapless runs are searched for gap minima, a touching counts only if its refined gap
converges below 1e-3 t, and it is a Weyl node only if its Berry charge is nonzero.

Referee: a 56^3 brute-force gap scan; each touching (gap < 1e-3 after refinement) is classified
by its Berry charge through spheres of two radii. If the radii disagree the case is inconclusive.
(An earlier referee used a fine Chern-number count; it flickers in exactly the hard cases.)
A proposal scores when the page's Weyl-node count differs from the referee's. Near-duplicates
count once.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from capture import cli  # noqa: E402
from kagome_core import WeylParams, weyl_berry_charge, weyl_detect_v2, weyl_gap_scan  # noqa: E402

TIMEOUT_S = 900
MAX_CASES = 12
BOUNDS = {"t": (0.2, 2.0), "breath": (-0.45, 0.45), "flux": (0.0, 360.0), "mA": (-3.0, 3.0),
          "tz": (-1.0, 1.0), "dz": (-1.5, 1.5)}
DEDUP = 0.03  # distance in canonical coordinates below which two failures count as one


def canonical(p):
    """Coordinates that can change node counts.

    tz shifts all three sublattices by the same 2*tz*cos(kz), so it never changes a gap and is
    dropped. Flux and 360-flux are mirror images (complex conjugation, k -> -k). Scaling t, mA, dz
    together only rescales the spectrum, so mA and dz are measured in units of t. An earlier
    version deduplicated on raw parameters, and an agent scored 12 by sweeping tz across one failure.
    """
    f = p.flux % 360.0
    f = min(f, 360.0 - f)
    return np.array([p.breath / 0.9, f / 180.0, p.mA / p.t / 6.0, p.dz / p.t / 3.0])


def check_case(c):
    if not isinstance(c, dict):
        return None, "case must be a dict"
    gap = c.get("gap", "12")
    if gap not in ("12", "23"):
        return None, "gap must be '12' or '23'"
    params = {k: v for k, v in c.items() if k != "gap"}
    for k, v in params.items():
        if k not in BOUNDS:
            return None, f"unknown parameter {k}"
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            return None, f"{k} must be finite"
        lo, hi = BOUNDS[k]
        if not lo <= v <= hi:
            return None, f"{k}={v} outside [{lo}, {hi}]"
    return (WeylParams(**params), gap), None


def score(result):
    if not isinstance(result, (list, tuple)) or not result:
        return False, "solution() must return a non-empty list of cases", {"combined_score": 0.0}
    cases = list(result)[:MAX_CASES]
    rows, failures = [], []
    for c in cases:
        parsed, err = check_case(c)
        if err:
            rows.append({"case": c, "status": "invalid", "reason": err})
            continue
        p, gap = parsed
        lower = 0 if gap == "12" else 1
        page = weyl_detect_v2(p, gap)["n_nodes"]
        charges = []
        for n in weyl_gap_scan(p, lower, 56):
            if n[3] < 1e-3:
                q1 = weyl_berry_charge(p, n[0], n[1], n[2], lower)
                q2 = weyl_berry_charge(p, n[0], n[1], n[2], lower, r=0.01, nt=16, nph=32)
                charges.append(q1 if q1 == q2 else None)
        ref_bf = sum(1 for q in charges if q)
        row = {"case": c, "page": page, "referee_nodes": ref_bf, "referee_charges": charges}
        if any(q is None for q in charges):
            row["status"] = "inconclusive (Berry charge depends on sphere radius)"
        elif page == ref_bf:
            row["status"] = "page correct"
        else:
            row["status"] = "page misses nodes" if page < ref_bf else "page reports false nodes"
            v = canonical(p)
            if all(np.linalg.norm(v - f["v"]) > DEDUP for f in failures if f["gap"] == gap):
                failures.append({"v": v, "gap": gap, "kind": row["status"]})
                row["counted"] = True
        rows.append(row)
    kinds = {f["kind"] for f in failures}
    s = len(failures) + (0.5 if len(kinds) == 2 else 0.0)
    public = {"cases": len(cases), "distinct_failures": len(failures),
              "misses": sum(f["kind"] == "page misses nodes" for f in failures),
              "false_nodes": sum(f["kind"] == "page reports false nodes" for f in failures),
              "inconclusive": sum(r.get("status", "").startswith("inconclusive") for r in rows)}
    return True, "evaluated", {"combined_score": float(s), "public": public, "private": {"rows": rows}}


if __name__ == "__main__":
    cli(score, TIMEOUT_S)
