"""Compare the Weyl page's node detector (run in Node from the HTML source) with a numpy referee.

Cases: the 40 random parameter sets of check_claims.py (W2, same seed) plus the cases the
AI-Scientist-v3 red team found. Referee per case: a 56^3 brute-force gap scan (touchings with
gap < 1e-3 after refinement), each touching classified by its Berry charge through spheres of two
radii (disagreement marks the case ambiguous). A fine Chern-number count is not used: it flickers
in exactly the hard cases. The page passes when its Weyl-node count equals the referee's.

Usage: python analysis/validate_weyl_detector.py   (needs node on PATH; about 20 min cold, referee results are cached)
Writes analysis/weyl_detector_validation.json.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from kagome_core import WeylParams, weyl_berry_charge, weyl_detect_like_js, weyl_gap_scan  # noqa: E402

cases = []
rng = np.random.default_rng(1)  # same stream as check_claims.py W2 (5 draws per trial)
for _ in range(40):
    cases.append({"breath": float(rng.uniform(-.3, .3)), "flux": float(rng.uniform(0, 360)), "mA": float(rng.uniform(-2, 2)),
                  "tz": float(rng.uniform(-.6, .6)), "dz": float(rng.uniform(-1.5, 1.5)), "source": "random"})
# found by the AI-Scientist-v3 red team (cycle 3): the gap-guard miss, the Chern-flicker false
# nodes, and two cases its referees disagreed on (both turn out to hold real Weyl nodes that a
# fine Chern-number count also misses)
cases += [
    {"breath": 0.08479690148362495, "flux": 306.94782185303643, "mA": 0.37176407241713605,
     "tz": -0.28788306271533215, "dz": 1.0196445630942264, "source": "agent: miss"},
    {"breath": 0.0, "flux": 1.0, "mA": -3.0, "tz": 0.0, "dz": 1.0, "source": "agent: false nodes"},
    {"breath": 0.0, "flux": 1.0, "mA": -3.0, "tz": 0.0, "dz": 0.8, "source": "agent: referee split"},
    {"breath": 0.0, "flux": 1.0, "mA": -2.5, "tz": 0.0, "dz": 0.3, "source": "agent: referee split"},
]

js_in = [{k: v for k, v in c.items() if k != "source"} | {"t": 1, "gapSel": "12"} for c in cases]
with tempfile.TemporaryDirectory() as td:
    fin, fout = Path(td) / "in.json", Path(td) / "out.json"
    fin.write_text(json.dumps(js_in))
    subprocess.run(["node", str(HERE / "run_page_weyl.js"), str(ROOT / "twisted-kagome-weyl-hyperlattice.html"), str(fin), str(fout)], check=True)
    page = json.loads(fout.read_text())

CACHE = HERE / "weyl_referee_cache.json"
cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}


def referee(p):
    """Touchings from a 56^3 scan, each classified by Berry charge at two sphere radii."""
    key = json.dumps([p.breath, p.flux, p.mA, p.tz, p.dz])
    if key not in cache:
        out = []
        for n in weyl_gap_scan(p, 0, 56):
            if n[3] >= 1e-3:
                continue
            q1 = weyl_berry_charge(p, n[0], n[1], n[2], 0)
            q2 = weyl_berry_charge(p, n[0], n[1], n[2], 0, r=0.01, nt=16, nph=32)
            out.append({"kz": n[2], "gap": n[3], "chi": q1 if q1 == q2 else None})
        cache[key] = out
        CACHE.write_text(json.dumps(cache))
    return cache[key]


rows = []
for c, pg in zip(cases, page):
    p = WeylParams(t=1, **{k: v for k, v in c.items() if k != "source"})
    ref = referee(p)
    ambiguous = any(t["chi"] is None for t in ref)
    ref_nodes = sum(1 for t in ref if t["chi"])
    ref_neutral = sum(1 for t in ref if t["chi"] == 0)
    old = weyl_detect_like_js(p, "12")["n_nodes"]
    nodes = len(pg["nodes"])
    rows.append({**c, "ref_nodes": ref_nodes, "ref_neutral": ref_neutral, "ambiguous": ambiguous, "old_page": old,
                 "page_nodes": nodes, "page_neutral": pg["neutral"], "flicker_rejected": pg["flicker"],
                 "extended": pg["extended"], "new_ok": nodes == ref_nodes, "old_ok": old == ref_nodes, "ms": pg["ms"]})
    r = rows[-1]
    print(f"{c['source']:<22} ref {ref_nodes} (+{ref_neutral} neutral{', AMBIGUOUS' if ambiguous else ''}) | old {old:>2} {'ok' if r['old_ok'] else 'X '} | "
          f"new {nodes:>2} {'ok' if r['new_ok'] else 'X '} | flicker rejected {pg['flicker']}", flush=True)

summary = {"cases": len(rows), "ambiguous": sum(r["ambiguous"] for r in rows),
           "old_correct": sum(r["old_ok"] for r in rows), "new_correct": sum(r["new_ok"] for r in rows),
           "mean_ms_new": float(np.mean([r["ms"] for r in rows]))}
print(summary)
(HERE / "weyl_detector_validation.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
