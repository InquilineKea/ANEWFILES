#!/bin/bash
# Run AI-Scientist-v3 natively (no Docker) on one kagome task with small budgets.
#
#   bash agsci_tasks/run_native.sh smoke              # one tiny agent call, checks auth
#   bash agsci_tasks/run_native.sh kagome_flatband    # full pipeline, reduced budgets
#   bash agsci_tasks/run_native.sh kagome_redteam
#
# Needs AGSCI_DIR (an AI-Scientist-v3 checkout after `uv sync`) and HEADLESS_DIR (headless-cli
# built at the revision AI-Scientist-v3 pins). Agents run with full tool permissions in the
# current machine, so only use this inside a disposable environment.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
AGSCI_DIR=${AGSCI_DIR:?set AGSCI_DIR to the AI-Scientist-v3 checkout}
HEADLESS_DIR=${HEADLESS_DIR:?set HEADLESS_DIR to the built headless-cli checkout}
PIN="git+https://github.com/RobertTLange/headless-cli.git#93cd9b06b85f848af1308c41e018991b33907c5e"
REAL_NPX=$(command -v npx)

# npm here cannot fetch git specs, so point the pinned spec at the local build of the same revision.
SHIM=$(mktemp -d)
cat > "$SHIM/npx" <<EOF
#!/bin/bash
if [ "\$1" = "-y" ] && [ "\$2" = "$PIN" ]; then shift 2; exec node "$HEADLESS_DIR/dist/cli.js" "\$@"; fi
exec "$REAL_NPX" "\$@"
EOF
cat > "$SHIM/headless" <<EOF
#!/bin/bash
exec node "$HEADLESS_DIR/dist/cli.js" "\$@"
EOF
chmod +x "$SHIM/npx" "$SHIM/headless"
export PATH="$SHIM:$PATH" HEADLESS_BIN="$SHIM/headless"

target=${1:?usage: run_native.sh smoke|kagome_flatband|kagome_redteam}
if [ "$target" = smoke ]; then
  exec headless claude --usage --allow yolo --prompt "Reply with exactly: pong"
fi

cd "$AGSCI_DIR"
exec .venv/bin/agsci run "$REPO/agsci_tasks/$target" \
  container.enabled=false container.env_file='' \
  investigation.execution.max_idea_cycles="${IDEA_CYCLES:-3}" \
  ablations.execution.num_replications=1 \
  ablations.execution.max_evals_per_variant=2 \
  posthoc.execution.max_selected_analyses=1 \
  ${PHASES:+execution.phases=$PHASES}
