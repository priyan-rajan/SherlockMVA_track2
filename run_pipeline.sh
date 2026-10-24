#!/usr/bin/env bash
# Reproducible entry point for the BUB1B N1002K unified pipeline.
# Stage 1 (structure/druggability) + Stage 2 (WT vs N1002K MD).
#
# Usage:
#   bash run_pipeline.sh                 # full run (default: 10 ns MD x2)
#   bash run_pipeline.sh -profile test   # fast smoke test (~minutes)
#   bash run_pipeline.sh -profile structure_only   # skip MD
#   bash run_pipeline.sh -resume         # resume an interrupted run
#   bash run_pipeline.sh --md_ns 5       # override any params.* this way
#
# Requires: nextflow (>=24), micromamba, and the two envs below (auto-created).
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
mkdir -p results

ensure_env () {
  local name="$1" file="$2"
  if ! micromamba env list 2>/dev/null | awk '{print $1}' | grep -qx "$name"; then
    echo "[run_pipeline] creating env '$name' from $file"
    micromamba create -y -f "$file" -n "$name"
  else
    echo "[run_pipeline] env '$name' present"
  fi
}

ensure_env bubdock environment_struct.yml
ensure_env bubmd   environment_md.yml

nextflow run main.nf -ansi-log false "$@" 2>&1 | tee -a results/run.log
