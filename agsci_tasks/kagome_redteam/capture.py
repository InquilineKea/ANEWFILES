"""Shared evaluator plumbing for the kagome tasks (copied into each task as harness.py).

Follows the AI-Scientist-v3 evaluator contract: evaluate.py --program_path P --results_dir D
writes D/metrics.json (with combined_score) and D/correct.json, plus stdout/stderr logs.
"""
import contextlib
import importlib.util
import io
import json
import multiprocessing as mp
import os
import sys
import traceback


class Tee(io.TextIOBase):
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()
        return len(s)


def _load(program_path, fn="solution"):
    d = os.path.dirname(os.path.abspath(program_path))
    if d not in sys.path:
        sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location("program_module", program_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["program_module"] = mod
    spec.loader.exec_module(mod)
    return getattr(mod, fn)


def _worker(path, q):
    try:
        q.put({"ok": True, "result": _load(path)()})
    except Exception as e:  # noqa: BLE001
        q.put({"ok": False, "error": str(e), "traceback": traceback.format_exc()})


def run_program(path, timeout):
    ctx = mp.get_context("fork" if "fork" in mp.get_all_start_methods() else "spawn")
    q = ctx.Queue(maxsize=1)
    pr = ctx.Process(target=_worker, args=(path, q))
    pr.start()
    pr.join(timeout)
    if pr.is_alive():
        pr.kill()
        pr.join(2)
        raise TimeoutError(f"program exceeded {timeout} s")
    payload = q.get(timeout=5)
    if not payload["ok"]:
        raise RuntimeError(payload["error"] + "\n" + payload["traceback"])
    return payload["result"]


def evaluate_main(program_path, results_dir, score_fn, timeout, verbose=True):
    """score_fn(result) -> (correct: bool, message: str, metrics: dict with combined_score)."""
    os.makedirs(results_dir, exist_ok=True)
    out_f = open(os.path.join(results_dir, "stdout.log"), "w")
    err_f = open(os.path.join(results_dir, "stderr.log"), "w")
    so = Tee(out_f, sys.__stdout__) if verbose else out_f
    se = Tee(err_f, sys.__stderr__) if verbose else err_f
    correct, msg, metrics = False, None, {"combined_score": 0.0}
    with contextlib.redirect_stdout(so), contextlib.redirect_stderr(se):
        try:
            print(f"Evaluating {program_path}")
            result = run_program(program_path, timeout)
            correct, msg, metrics = score_fn(result)
            print(f"correct={correct} combined_score={metrics.get('combined_score')}")
        except Exception as e:  # noqa: BLE001
            msg = f"Execution error: {e}"
            print(msg, file=sys.stderr)
            traceback.print_exc()
            metrics = {"combined_score": 0.0, "error": msg}
    with open(os.path.join(results_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2, default=float)
    with open(os.path.join(results_dir, "correct.json"), "w") as f:
        json.dump({"correct": bool(correct), "message": msg}, f, indent=2)
    out_f.close()
    err_f.close()
    return metrics


def cli(score_fn, timeout):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--program_path", default="program_0/src/program.py")
    ap.add_argument("--results_dir", default="results")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    evaluate_main(a.program_path, a.results_dir, score_fn, timeout, verbose=not a.quiet)
