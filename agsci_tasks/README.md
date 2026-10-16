# AI-Scientist-v3 tasks for the kagome pages

Two task directories in the format of [SakanaAI/AI-Scientist-v3](https://github.com/SakanaAI/AI-Scientist-v3)
(`agsci`). Both pass `agsci doctor` ("scaffold ready") and `agsci run --dry-run`.

| Task | Use | Score | Starter score |
|---|---|---|---|
| `kagome_flatband` | diversity: find different mechanisms that isolate a flat moiré band | log10(gap / width) of the best band isolated by ≥ 0.01 t | 1.24 (ratio 17, AA-confined) |
| `kagome_redteam` | evaluation: find where the Weyl page's node detector is wrong | distinct verified failures (+0.5 for both kinds) | 0 (targets the fixed detector) |

Evaluators re-solve whatever the program returns with their own sealed copy of
`analysis/kagome_core.py`, so a program cannot report a score it did not earn.
After editing the core, run `./sync_core.sh` to refresh the copies.

## Running

From an AI-Scientist-v3 checkout (`uv sync`, then build its Docker image as its README says):

```bash
agsci doctor /path/to/ANEWFILES/agsci_tasks/kagome_flatband
agsci run /path/to/ANEWFILES/agsci_tasks/kagome_flatband \
  investigation.execution.max_idea_cycles=4 \
  ablations.execution.num_replications=1 ablations.execution.max_evals_per_variant=2 \
  posthoc.execution.max_selected_analyses=1
```

It needs agent credentials in the checkout's `.env` (e.g. `CLAUDE_CODE_OAUTH_TOKEN` or
`ANTHROPIC_API_KEY`; Bedrock via `AWS_*`). Agents run with full tool permissions, so keep the
Docker default. Both tasks are CPU-only; a single evaluation takes about 1 s.

To evaluate a program by hand:

```bash
python kagome_flatband/evaluate.py --program_path kagome_flatband/program_0/src/program.py --results_dir /tmp/out
```
