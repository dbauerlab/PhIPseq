#!/bin/bash
# Run run_sample.sh for samples in samples.tsv, MAX_PARALLEL at a time
# (each sample uses ~6 threads at peak: bowtie2 -p4 + cutadapt -j2).
# macOS ships bash 3.2 (no `wait -n`), so concurrency is batched in fixed-size
# groups rather than polling a live job count.
#
# Usage: run_all.sh [tc_id ...]   -- restrict to specific samples (default: all)
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLES_TSV="$PROJECT_DIR/scripts/samples.tsv"
MAX_PARALLEL=1  # the network share (SMB mount) errors under concurrent read streams -- run sequentially

mkdir -p "$PROJECT_DIR/results"

declare -a WANT=("$@")

matches_want() {
  local tc=$1
  [ ${#WANT[@]} -eq 0 ] && return 0
  for w in "${WANT[@]}"; do
    [ "$w" = "$tc" ] && return 0
  done
  return 1
}

declare -a batch_tc=()
launch_batch() {
  [ ${#batch_tc[@]} -eq 0 ] && return
  wait "${batch_pids[@]}"
  batch_tc=()
  batch_pids=()
}

declare -a batch_pids=()
while IFS=$'\t' read -r tc_id library sample_prefix; do
  matches_want "$tc_id" || continue
  echo "=== launching $tc_id ($library, $sample_prefix) ==="
  bash "$PROJECT_DIR/scripts/run_sample.sh" "$tc_id" "$library" "$sample_prefix" \
    > "$PROJECT_DIR/results/${tc_id}.run.log" 2>&1 &
  batch_pids+=($!)
  batch_tc+=("$tc_id")
  if [ "${#batch_tc[@]}" -ge "$MAX_PARALLEL" ]; then
    launch_batch
  fi
done < <(tail -n +2 "$SAMPLES_TSV")

launch_batch
echo "=== all requested samples done ==="
