#!/bin/sh
# Copy the shared physics core and harness into each self-contained task directory.
# AI-Scientist-v3 snapshots a task directory, so tasks cannot import from outside it.
set -e
here=$(cd "$(dirname "$0")" && pwd)
for t in kagome_flatband kagome_redteam; do
  cp "$here/../analysis/kagome_core.py" "$here/$t/kagome_core.py"          # sealed copy used by the evaluator
  cp "$here/../analysis/kagome_core.py" "$here/$t/program_0/src/kagome_core.py"  # editable copy for the agent
  cp "$here/_harness.py" "$here/$t/capture.py"
done
echo synced
