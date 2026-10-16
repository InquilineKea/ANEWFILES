"""Evaluator: isolated flat band in a commensurate twisted kagome bilayer.

The program returns model parameters. The evaluator re-solves the supercell itself, so a
program cannot report a score it did not earn. combined_score = log10(gap / width) of the
best isolated band whose gap to both neighbours is at least MIN_GAP (in units of t).
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from capture import cli  # noqa: E402
from kagome_core import (MoireParams, isolated_band_scores, ldos_window, localisation_metrics,  # noqa: E402
                         moire_bands, moire_build, moire_grid, moire_path)

TIMEOUT_S = 600
MIN_GAP = 0.01
BOUNDS = {"m": (1, 4), "t": (1.0, 1.0), "breath": (-0.4, 0.4), "flux": (0.0, 360.0), "tp": (-1.0, 1.0),
          "xi": (0.08, 0.6), "rc": (0.2, 1.2), "V": (-1.0, 1.0), "sx": (0.0, 1.0), "sy": (0.0, 1.0)}


def validate(params):
    if not isinstance(params, dict):
        return False, "solution() must return a dict of model parameters"
    unknown = set(params) - set(BOUNDS)
    if unknown:
        return False, f"unknown parameters: {sorted(unknown)}"
    for k, v in params.items():
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            return False, f"{k} must be a finite number"
        lo, hi = BOUNDS[k]
        if not lo - 1e-12 <= v <= hi + 1e-12:
            return False, f"{k}={v} outside [{lo}, {hi}]"
    if "m" in params and int(params["m"]) != params["m"]:
        return False, "m must be an integer"
    return True, "parameters valid"


def score(result):
    ok, msg = validate(result)
    if not ok:
        return False, msg, {"combined_score": 0.0, "error": msg}
    p = MoireParams(**{**result, "m": int(result.get("m", 3))})
    sc = moire_build(p)
    ks = np.concatenate([moire_path(sc, 12), moire_grid(sc, 6)])
    E = moire_bands(sc, ks)
    bands = [b for b in isolated_band_scores(E) if b["gap"] >= MIN_GAP]
    if not bands:
        return True, "no band is isolated by at least MIN_GAP", {"combined_score": -3.0, "public": {"isolated_bands": 0}}
    best = max(bands, key=lambda b: b["ratio"])
    G = sc.high_symmetry()["G"]
    lo, hi = E[:, best["band"]].min(), E[:, best["band"]].max()
    wt, nst = ldos_window(sc, G, (lo + hi) / 2, (hi - lo) / 2 + 1e-9)
    loc = localisation_metrics(sc, wt)
    m0 = sc.n - 2 * sc.N
    public = {
        "theta_deg": math.degrees(sc.theta), "n_sites": sc.n, "isolated_bands": len(bands),
        "best_band_index": best["band"], "best_band_in_top_manifold": best["band"] >= m0,
        "best_band_energy": best["E"], "best_band_width": best["width"], "best_band_gap": best["gap"],
        "flatness_ratio": best["ratio"], "aa_enhancement": loc["aa_enhancement"],
        "participation": loc["participation"], "layer1_share": loc["layer1_share"],
    }
    return True, "evaluated", {"combined_score": float(math.log10(best["ratio"])), "public": public,
                               "private": {"params": result, "k_points": len(ks)}}


if __name__ == "__main__":
    cli(score, TIMEOUT_S)
