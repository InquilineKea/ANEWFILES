"""Compare the Weyl page's node detector (run in Node from the HTML source) with a numpy referee.

Cases: the 40 random parameter sets of check_claims.py (W2, same seed) plus the cases the
AI-Scientist-v3 red team found. Referee per case: a 56^3 brute-force gap scan (touchings with
gap < 1e-3 after refinement) and a fine Chern detector (30x30 grid, 97 kz planes). The page
reports Weyl nodes (nonzero Berry charge) and zero-charge touchings separately; nodes are scored
against the fine Chern count, nodes + touchings against the scan.

Usage: python analysis/validate_weyl_detector.py   (needs node on PATH; about 15 min)
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
from kagome_core import WeylParams, weyl_detect_like_js, weyl_gap_scan  # noqa: E402

cases = []
rng = np.random.default_rng(1)  # same stream as check_claims.py W2 (5 draws per trial)
for _ in range(40):
    cases.append({"breath": float(rng.uniform(-.3, .3)), "flux": float(rng.uniform(0, 360)), "mA": float(rng.uniform(-2, 2)),
                  "tz": float(rng.uniform(-.6, .6)), "dz": float(rng.uniform(-1.5, 1.5)), "source": "random"})
# found by the AI-Scientist-v3 red team (cycle 3): the gap-guard miss, the Chern-flicker false
# nodes, and two cases with charge-free touchings
cases += [
    {"breath": 0.08479690148362495, "flux": 306.94782185303643, "mA": 0.37176407241713605,
     "tz": -0.28788306271533215, "dz": 1.0196445630942264, "source": "agent: miss"},
    {"breath": 0.0, "flux": 1.0, "mA": -3.0, "tz": 0.0, "dz": 1.0, "source": "agent: false nodes"},
    {"breath": 0.0, "flux": 1.0, "mA": -3.0, "tz": 0.0, "dz": 0.8, "source": "agent: inconclusive"},
    {"breath": 0.0, "flux": 1.0, "mA": -2.5, "tz": 0.0, "dz": 0.3, "source": "agent: inconclusive"},
]

js_in = [{k: v for k, v in c.items() if k != "source"} | {"t": 1, "gapSel": "12"} for c in cases]
with tempfile.TemporaryDirectory() as td:
    fin, fout = Path(td) / "in.json", Path(td) / "out.json"
    fin.write_text(json.dumps(js_in))
    subprocess.run(["node", str(HERE / "run_page_weyl.js"), str(ROOT / "twisted-kagome-weyl-hyperlattice.html"), str(fin), str(fout)], check=True)
    page = json.loads(fout.read_text())

rows = []
for c, pg in zip(cases, page):
    p = WeylParams(t=1, **{k: v for k, v in c.items() if k != "source"})
    scan = len([n for n in weyl_gap_scan(p, 0, 56) if n[3] < 1e-3])
    chern = weyl_detect_like_js(p, "12", nk=30, nz=96)["n_nodes"]
    old = weyl_detect_like_js(p, "12")["n_nodes"]
    nodes, touch = len(pg["nodes"]), pg["neutral"]
    ok_nodes = nodes == chern
    ok_total = nodes + touch == scan
    rows.append({**c, "referee_scan": scan, "referee_chern": chern, "old_page": old, "page_nodes": nodes,
                 "page_neutral": touch, "flicker_rejected": pg["flicker"], "extended": pg["extended"],
                 "ok": bool(ok_nodes and ok_total), "ms": pg["ms"]})
    print(f"{c['source']:<22} old {old:>2} | new nodes {nodes:>2} + touchings {touch} | referee chern {chern:>2}, scan {scan:>2} | "
          f"{'ok' if rows[-1]['ok'] else 'MISMATCH'}", flush=True)

summary = {
    "cases": len(rows),
    "new_correct": sum(r["ok"] for r in rows),
    "old_correct": sum(r["old_page"] == r["referee_scan"] == r["referee_chern"] or r["old_page"] == r["referee_chern"] for r in rows),
    "mean_ms": float(np.mean([r["ms"] for r in rows])),
}
print(summary)
(HERE / "weyl_detector_validation.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
