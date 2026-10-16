"""Evaluator: red-team the Weyl-node detector of the Twisted Kagome-Weyl Hyperlattice page.

The page counts Weyl nodes between two bands from jumps in the Chern number C(kz), using an
18x18 Chern grid, 49 kz planes, and skipping planes whose gap is below 0.01 t. The program
proposes parameter sets; each one is checked by a referee that must agree with itself at two
resolutions (a 30x30 / 97-plane detector and a 56^3 brute-force gap scan). A proposal scores
when the page's count differs from an agreed referee count. Near-duplicates count once.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from capture import cli  # noqa: E402
from kagome_core import WeylParams, weyl_detect_like_js, weyl_gap_scan  # noqa: E402

TIMEOUT_S = 900
MAX_CASES = 12
BOUNDS = {"t": (0.2, 2.0), "breath": (-0.45, 0.45), "flux": (0.0, 360.0), "mA": (-3.0, 3.0),
          "tz": (-1.0, 1.0), "dz": (-1.5, 1.5)}
SCALE = np.array([1.8, 0.9, 360.0, 6.0, 2.0, 3.0])
DEDUP = 0.03  # normalised parameter distance below which two failures count as one


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
        page = weyl_detect_like_js(p, gap)["n_nodes"]
        ref_det = weyl_detect_like_js(p, gap, nk=30, nz=96)["n_nodes"]
        ref_bf = len([n for n in weyl_gap_scan(p, lower, 56) if n[3] < 1e-3])
        row = {"case": c, "page": page, "referee_detector": ref_det, "referee_scan": ref_bf}
        if ref_det != ref_bf:
            row["status"] = "inconclusive (referees disagree)"
        elif page == ref_bf:
            row["status"] = "page correct"
        else:
            row["status"] = "page misses nodes" if page < ref_bf else "page reports false nodes"
            v = np.array([p.t, p.breath, p.flux, p.mA, p.tz, p.dz]) / SCALE
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
