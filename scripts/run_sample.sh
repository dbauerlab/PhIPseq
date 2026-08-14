#!/bin/bash
# Trim + align + count one sample against the combined bowtie2 index, in
# paired-end mode.
#
# Both sequencing lanes are merged into a single stream per mate *before*
# trimming/alignment (rather than processing each lane separately and
# combining counts afterward). cutadapt trims R1/R2 together in one paired
# invocation, then bowtie2 aligns each read pair jointly against the index
# (rather than trimming/aligning R1 and R2 independently as single-end
# reads, which is what this script did before the SLURM/HPC migration).
#
# Usage: run_sample.sh <tc_id> <library: CoV|Vir3> <sample_prefix> [subsample_reads]
#   subsample_reads: if >0, only process the first N read pairs (taken from
#     the merged lanes) -- for a quick sanity-check run before committing to
#     the full data.
set -euo pipefail

TC_ID=$1
LIBRARY=$2
SAMPLE_PREFIX=$3
SUBSAMPLE=${4:-0}

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDEX="$PROJECT_DIR/reference/bowtie2_index/combined"
OUT_DIR="$PROJECT_DIR/results/$TC_ID"

# A fresh, unique temp directory per invocation -- honors $TMPDIR if the
# scheduler sets one (SLURM typically provisions node-local scratch there),
# falls back to /tmp otherwise. Replaces the old fixed
# /tmp/phipseq_tmp/$TC_ID path, which risked collisions across concurrent
# array tasks/reruns. Cleaned up on exit even if the script fails partway
# through.
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

mkdir -p "$OUT_DIR"

# --- Load config ------------------------------------------------------------
# Paths/tunables/adapters all live in reference/*.json rather than being
# hardcoded here, so redeploying to a new machine only means editing those
# files. This one python3 call loads all three and prints shell variable
# assignments for `eval` to pick up -- python3 only needs the stdlib `json`
# module here, so it works whether or not the "phipseq" conda env is active.
eval "$(python3 - "$PROJECT_DIR" "$LIBRARY" <<'PY'
import json
import shlex
import sys

project_dir, library = sys.argv[1], sys.argv[2]
paths = json.load(open(f"{project_dir}/reference/paths.json"))
params = json.load(open(f"{project_dir}/reference/run_params.json"))
adapters = json.load(open(f"{project_dir}/reference/adapters.json"))[library]

print(f"FASTQ_DIR={shlex.quote(paths['fastq_dir'])}")
print(f"BOWTIE_THREADS={shlex.quote(str(params['bowtie_threads']))}")
print(f"CUTADAPT_THREADS={shlex.quote(str(params['cutadapt_threads']))}")
print("LANES=(" + " ".join(shlex.quote(lane) for lane in params["lanes"]) + ")")
print(f"R1_ADAPTER={shlex.quote(adapters['R1'])}")
print(f"R2_ADAPTER={shlex.quote(adapters['R2'])}")
PY
)"

# --- Step 1: merge lanes, per mate ------------------------------------------
# Concatenate every lane's raw FASTQ for one mate into a single file before
# any trimming/alignment happens. Gzip-concatenation is valid (a
# concatenated .gz is a valid multi-member gzip stream), so the full-data
# path just `cat`s the raw .gz lane files together -- no decompression
# needed here; cutadapt reads gzip input natively.
merge_mate() {
  local mate=$1 out=$2
  local lane_files=()
  for lane in "${LANES[@]}"; do
    lane_files+=("$FASTQ_DIR/${SAMPLE_PREFIX}_${lane}_${mate}_001.fastq.gz")
  done
  if [ "$SUBSAMPLE" -gt 0 ]; then
    # Decompress, merge, and truncate to the first SUBSAMPLE read pairs
    # (same line count taken from R1 and R2 keeps the two mates in sync),
    # then recompress so downstream steps only ever see .fastq.gz input.
    zcat "${lane_files[@]}" | head -n $((SUBSAMPLE * 4)) | gzip > "$out"
  else
    cat "${lane_files[@]}" > "$out"
  fi
}

R1_MERGED="$TMP_DIR/R1_merged.fastq.gz"
R2_MERGED="$TMP_DIR/R2_merged.fastq.gz"
echo "[$TC_ID] merging lanes (${LANES[*]}) per mate..." >&2
merge_mate R1 "$R1_MERGED"
merge_mate R2 "$R2_MERGED"

# --- Step 2: paired-mode trimming -------------------------------------------
# Trim R1 and R2 together in one cutadapt invocation so mate pairing is
# preserved, instead of two independent single-mate calls. cutadapt's
# default paired behavior discards a whole pair if either mate fails
# --discard-untrimmed (--pair-filter=any) -- accepted as-is.
#
# Output goes to real temp files (-o/-p), not stdout, because paired-mode
# bowtie2 needs two complete, separate inputs (-1/-2) rather than a single
# stdin stream -- a plain shell pipe can't hand off two synchronized FASTQ
# streams at once.
R1_TRIMMED="$TMP_DIR/R1_trimmed.fastq.gz"
R2_TRIMMED="$TMP_DIR/R2_trimmed.fastq.gz"
echo "[$TC_ID] trimming (paired mode)..." >&2
cutadapt -j "$CUTADAPT_THREADS" \
  -g "$R1_ADAPTER" -G "$R2_ADAPTER" \
  --discard-untrimmed \
  -o "$R1_TRIMMED" -p "$R2_TRIMMED" \
  "$R1_MERGED" "$R2_MERGED" \
  > "$OUT_DIR/cutadapt.log" 2>&1

# --- Step 3: paired-end alignment + counting --------------------------------
# Align mate pairs jointly (-1/-2), rather than each mate independently
# against the index. Default --fr orientation (mate1 forward / mate2
# reverse) is correct for this library prep, so no orientation override is
# needed here (the old single-end --norc/--nofw hack is gone).
#
# A paired-mode BAM has one alignment record per *mate*, so filtering on
# just "-F 4" (mapped) would double-count each fragment. "-f 64" (0x40,
# first-in-pair) keeps exactly one record per aligned pair, so the
# cut -f3 | sort | uniq -c counting step below still produces one count per
# fragment, same shape as before.
echo "[$TC_ID] aligning (paired mode)..." >&2
bowtie2 -p "$BOWTIE_THREADS" -x "$INDEX" \
  -1 "$R1_TRIMMED" -2 "$R2_TRIMMED" \
  --no-unal \
  2> "$OUT_DIR/bowtie2.log" \
  | samtools view -F 4 -f 64 - \
  | cut -f3 \
  | sort | uniq -c | awk '{print $2"\t"$1}' > "$OUT_DIR/counts.tsv"

echo "[$TC_ID] done." >&2
