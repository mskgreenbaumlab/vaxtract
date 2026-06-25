#!/usr/bin/env bash
# Pre-flight: prove a COMPUTE node can (a) reach the Claude API and (b) use the subscription login,
# from INSIDE the container, before launching a batch. Run it under srun on a compute node, e.g.:
#
#   srun -A your-slurm-account -p componc_cpu -t 00:10:00 --mem 4G \
#        bash scale/probe_egress.sh /path/to/containers/antvac_extractor.sif $HOME/.claude
#
# Exit 0 = egress + auth OK (safe to run the batch). Non-zero = fix that before scaling.
set -euo pipefail

SIF="${1:?usage: probe_egress.sh <sif> <creds_dir>}"
CREDS="${2:?usage: probe_egress.sh <sif> <creds_dir>}"
MOUNT="/root/.claude"

echo "== 1. raw HTTPS egress to api.anthropic.com =="
singularity exec --bind /path/to/data "$SIF" \
    bash -lc 'curl -sS -o /dev/null -w "  HTTP %{http_code} in %{time_total}s\n" https://api.anthropic.com/v1/ || { echo "  NO EGRESS"; exit 2; }'

echo "== 2. claude CLI sees the subscription login + completes a tiny prompt =="
# Mirror the batch's auth mode: extraction_agent.py --subscription drops ANTHROPIC_API_KEY before
# spawning claude, so the CLI uses the plan quota. The host commonly exports a (possibly depleted)
# ANTHROPIC_API_KEY via .bashrc, which singularity leaks into the container; if we DON'T unset it
# here the probe bills against API credit and reports "Credit balance is too low" — a false negative
# that does not reflect the real subscription run. So unset it for this check.
singularity exec --bind "/path/to/data,${CREDS}:${MOUNT}" --env "CLAUDE_CONFIG_DIR=${MOUNT}" "$SIF" \
    bash -lc 'unset ANTHROPIC_API_KEY; claude --version && echo "say OK" | claude -p "reply with exactly: OK"'

echo "== probe passed: egress + subscription auth work from this node =="
