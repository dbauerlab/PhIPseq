#!/bin/bash
# Trim + align + count one sample (both lanes, both mates) against the
# combined bowtie2 index. Streams cutadapt -> bowtie2 -> samtools with no
# intermediate files, so nothing large ever touches the (Dropbox-synced)
# project directory.
#
# Usage: run_sample.sh <tc_id> <library: CoV|Vir3> <sample_prefix> [subsample_reads]
#   subsample_reads: if >0, only process the first N reads per lane per mate
#     (for a quick sanity-check run before committing to the full data).
set -euo pipefail

TC_ID=$1
LIBRARY=$2
SAMPLE_PREFIX=$3
SUBSAMPLE=${4:-0}

FASTQ_DIR="/Volumes/lab-bauerd/data/STPs/genomics-stp/inputs/tiphaine.cayol/PM26049/20260810_LH00442_0282_A23G2Y5LT3/fastq"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDEX="$PROJECT_DIR/reference/bowtie2_index/combined"
ADAPTERS_JSON="$PROJECT_DIR/reference/adapters.json"
OUT_DIR="$PROJECT_DIR/results/$TC_ID"
TMP_DIR="/tmp/phipseq_tmp/$TC_ID"  # local, not Dropbox-synced -- intermediate hit lists live here
LANES=(L007 L008)
BOWTIE_THREADS=4
CUTADAPT_THREADS=2

mkdir -p "$OUT_DIR" "$TMP_DIR"

R1_ADAPTER=$("$PROJECT_DIR/.venv/bin/python" -c "import json;print(json.load(open('$ADAPTERS_JSON'))['$LIBRARY']['R1'])")
R2_ADAPTER=$("$PROJECT_DIR/.venv/bin/python" -c "import json;print(json.load(open('$ADAPTERS_JSON'))['$LIBRARY']['R2'])")

read_fastq() {
  # `cat | gzcat` rather than `gzcat` directly on the network file: some files
  # on this SMB mount reproducibly fail close() with EBADF even though every
  # byte reads back correctly (confirmed: full read count matches file size
  # exactly). `cat` doesn't check close()'s return value so it's unaffected;
  # piping means gzcat decompresses from a plain pipe, never touching the
  # problematic file handle itself.
  # `|| true`: head closing early on a subsample sends the pipe SIGPIPE,
  # which pipefail would otherwise propagate as a spurious failure.
  if [ "$SUBSAMPLE" -gt 0 ]; then
    { cat "$1" 2>/dev/null | gzcat 2>/dev/null | head -n $((SUBSAMPLE * 4)); } || true
  else
    cat "$1" | gzcat
  fi
}

process_mate() {
  local mate=$1 adapter=$2 strand_flag=$3
  local lane_hits=()
  for lane in "${LANES[@]}"; do
    fq="$FASTQ_DIR/${SAMPLE_PREFIX}_${lane}_${mate}_001.fastq.gz"
    local hits="$TMP_DIR/${mate}_${lane}_hits.tmp"
    lane_hits+=("$hits")
    local ok=false
    for attempt in 1 2 3 4 5; do
      echo "[$TC_ID] $mate $lane: trimming + aligning (attempt $attempt)..." >&2
      # `> "$hits"` (not `>>`) inside the loop: a retry must overwrite, not
      # append, or a partial failed attempt would double-count on success.
      if read_fastq "$fq" \
           | cutadapt -j "$CUTADAPT_THREADS" -g "$adapter" --discard-untrimmed -o - - \
               2> "$OUT_DIR/cutadapt_${mate}_${lane}.log" \
           | bowtie2 -p "$BOWTIE_THREADS" -x "$INDEX" "$strand_flag" --no-unal -U - \
               2> "$OUT_DIR/bowtie2_${mate}_${lane}.log" \
           | samtools view -F 4 - | cut -f3 > "$hits"
      then
        ok=true
        break
      fi
      echo "[$TC_ID] $mate $lane: attempt $attempt failed (likely a transient SMB read error), retrying in 20s..." >&2
      sleep 20
    done
    if [ "$ok" != true ]; then
      echo "[$TC_ID] $mate $lane: FAILED after 5 attempts, giving up." >&2
      exit 1
    fi
  done
  cat "${lane_hits[@]}" | sort | uniq -c | awk '{print $2"\t"$1}' > "$OUT_DIR/${mate}_counts.tsv"
  rm -f "${lane_hits[@]}"
}

process_mate R1 "$R1_ADAPTER" --norc
process_mate R2 "$R2_ADAPTER" --nofw

echo "[$TC_ID] done." >&2
