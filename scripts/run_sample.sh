#!/bin/bash
# (The line above is the "shebang": it tells the OS to run this file with
# /bin/bash when it is executed directly, e.g. `./run_sample.sh ...`.)
#
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
#
# Pipeline outputs (all in <run_dir>/results/<tc_id>/):
#   cutadapt.log -- cutadapt's trimming report (read in by summarize.py for QC)
#   bowtie2.log  -- bowtie2's alignment-rate summary (read in by summarize.py)
#   counts.tsv   -- two tab-separated columns, no header: <ref_id> <count>,
#                   one row per peptide that at least one read pair hit.
#
# Strict mode, roughly "make bash behave more like Python with exceptions":
#   -e           exit immediately if any command fails (non-zero exit status)
#                instead of carrying on regardless (bash's default).
#   -u           treat use of an unset variable as an error (like a NameError)
#                rather than silently substituting an empty string.
#   -o pipefail  a pipeline `a | b | c` fails if *any* stage fails, not just
#                the last one (by default only `c`'s exit status counts).
set -euo pipefail

# Positional arguments: $1, $2, ... are the command-line arguments (like
# sys.argv[1], sys.argv[2], ...). With `set -u` above, a missing required
# argument makes the script stop with an "unbound variable" error.
TC_ID=$1
LIBRARY=$2
SAMPLE_PREFIX=$3
# ${4:-0} means "the 4th argument if given and non-empty, otherwise 0" --
# bash's way of giving an optional argument a default value.
SUBSAMPLE=${4:-0}

# Absolute path to the repo checkout (the parent of this scripts/ folder),
# Python equivalent: Path(__file__).resolve().parent.parent. Read inside-out:
#   ${BASH_SOURCE[0]}   path of this script file as it was invoked
#   dirname ...         strip the filename -> the scripts/ directory
#   $( ... )            "command substitution": run the command and paste its
#                       printed output in place (like subprocess output)
#   cd .../.. && pwd    move to the parent directory and print its absolute
#                       path; this runs in a subshell, so it does not change
#                       the current directory of the script itself.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# RUN_DIR (where results/ and the generated bowtie2 index live) is resolved
# below from reference/paths.json, once config is loaded -- it's a separate
# directory from PROJECT_DIR (the repo checkout) so the repo stays free of
# run output.

# A fresh, unique temp directory per invocation -- honors $TMPDIR if the
# scheduler sets one (SLURM typically provisions node-local scratch there),
# falls back to /tmp otherwise. Replaces the old fixed
# /tmp/phipseq_tmp/$TC_ID path, which risked collisions across concurrent
# array tasks/reruns. Cleaned up on exit even if the script fails partway
# through.
# `mktemp -d` creates a new uniquely named empty directory and prints its
# path (like Python's tempfile.mkdtemp()).
TMP_DIR="$(mktemp -d)"
# `trap '<command>' EXIT` registers a command to run whenever the script
# exits, for any reason (success, error under `set -e`, Ctrl-C) -- similar to
# a try/finally or atexit handler. Here it deletes the temp directory and
# everything in it. Single quotes delay expansion of $TMP_DIR until the trap
# actually fires.
trap 'rm -rf "$TMP_DIR"' EXIT

# --- Load config ------------------------------------------------------------
# Paths/tunables/adapters all live in reference/*.json (in the repo
# checkout) rather than being hardcoded here, so redeploying to a new
# machine only means editing those files. This one python3 call loads all
# three and prints shell variable assignments for `eval` to pick up --
# python3 only needs the stdlib `json` module here, so it works whether or
# not the "phipseq" conda env is active.
#
# How the next block works, piece by piece:
#   python3 - "$PROJECT_DIR" "$LIBRARY"
#       runs Python reading its program from stdin (the lone `-`); the two
#       extra words become sys.argv[1] and sys.argv[2] inside Python.
#   <<'PY' ... PY
#       a "heredoc": every line up to the closing `PY` line is fed to
#       python3's stdin as the program text. Quoting the marker ('PY') means
#       bash does NOT expand $variables or $(...) inside it, so the Python
#       code (f-strings with {braces} etc.) is passed through verbatim.
#   "$( ... )"
#       captures everything the Python program prints, e.g. lines like
#       FASTQ_DIR=/some/path and LANES=(L006) (values get quotes added
#       only when they contain spaces or special characters).
#   eval "..."
#       runs that captured text as bash code, so those lines become real
#       shell variables in this script. shlex.quote() on the Python side
#       guarantees each value is safely quoted for bash first.
eval "$(python3 - "$PROJECT_DIR" "$LIBRARY" <<'PY'
# (Python code from here to the closing PY marker.)
import json
import shlex
import sys

# The two arguments passed after "python3 -" above.
project_dir, library = sys.argv[1], sys.argv[2]
# External file locations + run_dir (where all output is written).
paths = json.load(open(f"{project_dir}/reference/paths.json"))
# Lanes to merge and thread counts for cutadapt/bowtie2.
params = json.load(open(f"{project_dir}/reference/run_params.json"))
# R1/R2 adapter sequences for just the library of this sample ("CoV" or "Vir3");
# an unknown library name raises KeyError. (eval itself still "succeeds" on
# the empty output, so the script then stops at the first use of an unset
# variable such as $RUN_DIR, thanks to set -u.)
adapters = json.load(open(f"{project_dir}/reference/adapters.json"))[library]

# Each print() emits one line of bash code, of the form NAME=value.
# Directory holding the raw *.fastq.gz files.
print(f"FASTQ_DIR={shlex.quote(paths['fastq_dir'])}")
# Root directory for generated reference data and results/.
print(f"RUN_DIR={shlex.quote(paths['run_dir'])}")
# Threads for bowtie2 (-p) and cutadapt (-j); numbers converted to str for quoting.
print(f"BOWTIE_THREADS={shlex.quote(str(params['bowtie_threads']))}")
print(f"CUTADAPT_THREADS={shlex.quote(str(params['cutadapt_threads']))}")
# A bash *array* literal, e.g. LANES=(L006) or LANES=(L001 L002) -- the
# bash equivalent of a Python list.
print("LANES=(" + " ".join(shlex.quote(lane) for lane in params["lanes"]) + ")")
# The constant 5-prime sequence to trim off the start of mate 1 / mate 2.
print(f"R1_ADAPTER={shlex.quote(adapters['R1'])}")
print(f"R2_ADAPTER={shlex.quote(adapters['R2'])}")
PY
)"

# RUN_DIR (separate from PROJECT_DIR, the repo checkout) is where generated
# reference data and results/ live -- see reference/paths.json.
# INDEX is the bowtie2 index "basename": bowtie2 looks for files named
# combined.1.bt2, combined.2.bt2, ... in that folder (built by
# build_index.sbatch).
INDEX="$RUN_DIR/reference/bowtie2_index/combined"
# Per-sample output folder for counts.tsv and the two log files.
OUT_DIR="$RUN_DIR/results/$TC_ID"
# mkdir -p: create the directory and any missing parents; no error if it
# already exists (like os.makedirs(..., exist_ok=True)).
mkdir -p "$OUT_DIR"

# --- Step 1: merge lanes, per mate ------------------------------------------
# Concatenate every lane's raw FASTQ for one mate into a single file before
# any trimming/alignment happens. Gzip-concatenation is valid (a
# concatenated .gz is a valid multi-member gzip stream), so the full-data
# path just `cat`s the raw .gz lane files together -- no decompression
# needed here; cutadapt reads gzip input natively.
#
# A bash function: called as `merge_mate R1 /path/out.fastq.gz`; inside it
# the arguments are $1 and $2 (functions have their own positional args).
merge_mate() {
  # `local` makes these variables private to the function (by default every
  # bash variable is global). mate = "R1" or "R2"; out = output file path.
  local mate=$1 out=$2
  # An empty local array to collect this mate's per-lane file paths.
  local lane_files=()
  # Loop over every lane. "${LANES[@]}" expands to each array element as a
  # separate, safely quoted word (like iterating over a Python list).
  for lane in "${LANES[@]}"; do
    # += appends an element to the array. Builds the Illumina-style name,
    # e.g. <FASTQ_DIR>/CAY9116A13_S4_L006_R1_001.fastq.gz. ${VAR} braces just
    # mark where a variable name ends next to other text such as "_".
    lane_files+=("$FASTQ_DIR/${SAMPLE_PREFIX}_${lane}_${mate}_001.fastq.gz")
  done
  # [ a -gt b ] is bash's integer "a > b" test.
  if [ "$SUBSAMPLE" -gt 0 ]; then
    # Decompress, merge, and truncate to the first SUBSAMPLE read pairs
    # (same line count taken from R1 and R2 keeps the two mates in sync),
    # then recompress so downstream steps only ever see .fastq.gz input.
    # Pipeline (each `|` feeds one command's stdout into the next's stdin):
    #   zcat <files...>   decompress all lane files, concatenated, to stdout
    #   head -n N         keep only the first N lines; $((...)) is bash
    #                     integer arithmetic -- a FASTQ record is 4 lines, so
    #                     SUBSAMPLE*4 lines = SUBSAMPLE reads
    #   gzip > "$out"     recompress and write (`>` redirects stdout to a
    #                     file, overwriting it)
    zcat "${lane_files[@]}" | head -n $((SUBSAMPLE * 4)) | gzip > "$out"
  else
    # Full data: byte-level concatenation of the gzipped lane files into one
    # gzip file (no decompression needed, see block comment above).
    cat "${lane_files[@]}" > "$out"
  fi
}

# Temp-dir paths for each mate's merged (but not yet trimmed) reads.
R1_MERGED="$TMP_DIR/R1_merged.fastq.gz"
R2_MERGED="$TMP_DIR/R2_merged.fastq.gz"
# Progress message. ${LANES[*]} joins the array into one space-separated
# string for display; `>&2` sends the message to stderr (so it lands in the
# SLURM job log, not mixed into any data stream on stdout).
echo "[$TC_ID] merging lanes (${LANES[*]}) per mate..." >&2
# Run the function above once per mate.
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
# Temp-dir paths for each mate's trimmed reads.
R1_TRIMMED="$TMP_DIR/R1_trimmed.fastq.gz"
R2_TRIMMED="$TMP_DIR/R2_trimmed.fastq.gz"
echo "[$TC_ID] trimming (paired mode)..." >&2
# The cutadapt command below spans several lines: a trailing `\` continues
# the same command on the next line. Its options:
#   -j N                  use N CPU cores (CUTADAPT_THREADS from run_params.json)
#   -g SEQ                a 5' ("front") adapter to find in mate 1 (R1); the
#                         adapter and everything before it is removed, leaving
#                         the peptide-encoding insert
#   -G SEQ                same, for mate 2 (R2) -- uppercase options act on R2
#   --discard-untrimmed   drop reads where the adapter was NOT found (in paired
#                         mode the whole pair is dropped, see above)
#   -o FILE / -p FILE     output files for trimmed R1 / trimmed R2 (gzip
#                         output chosen automatically from the .gz suffix)
#   "$R1_MERGED" "$R2_MERGED"  the two input files (mate 1, mate 2)
#   > cutadapt.log 2>&1   send stdout (the trimming report) to cutadapt.log,
#                         and `2>&1` sends stderr to the same place, so
#                         warnings/errors end up in that log too.
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
# (Strictly, bowtie2 writes uncompressed SAM text to stdout, not BAM;
# samtools view reads either.)
echo "[$TC_ID] aligning (paired mode)..." >&2
# One long pipeline, continued over several lines with trailing `\`:
#   bowtie2                 the aligner; writes alignments (SAM) to stdout
#     -p N                  use N threads (BOWTIE_THREADS from run_params.json)
#     -x INDEX              index basename to align against (see INDEX above)
#     -1 FILE -2 FILE       trimmed mate-1 and mate-2 reads (paired-end mode)
#     --no-unal             don't output records for reads that aligned nowhere
#     2> bowtie2.log        stderr only (bowtie2's alignment-rate summary) goes
#                           to bowtie2.log; stdout keeps flowing down the pipe
#   samtools view - ...     read SAM from stdin (`-`), print matching records
#                           as SAM text lines (no header, as -h isn't given)
#     -F 4                  exclude records with flag 4 set (= "read unmapped")
#     -f 64                 require flag 64 set (= "first mate of the pair")
#   cut -f3                 keep only tab-separated column 3 of each SAM line,
#                           the reference name = the peptide ref_id hit
#   sort                    sort the ref_ids so identical ones are adjacent
#                           (uniq only collapses *adjacent* duplicates)
#   uniq -c                 collapse duplicates, prefixing each with its count:
#                           "   123 ref_id"
#   awk '{print $2"\t"$1}'  swap the two whitespace-separated fields and join
#                           them with a tab: "ref_id<TAB>123" ($1/$2 here are
#                           awk's field numbers, not bash variables -- the
#                           single quotes stop bash from touching them)
#   > counts.tsv            write the final table to the sample's counts.tsv
# Because of `set -o pipefail`, a failure in any stage fails the script.
bowtie2 -p "$BOWTIE_THREADS" -x "$INDEX" \
  -1 "$R1_TRIMMED" -2 "$R2_TRIMMED" \
  --no-unal \
  2> "$OUT_DIR/bowtie2.log" \
  | samtools view -F 4 -f 64 - \
  | cut -f3 \
  | sort | uniq -c | awk '{print $2"\t"$1}' > "$OUT_DIR/counts.tsv"

# Final progress message (to stderr). The EXIT trap above then removes
# TMP_DIR, including the merged and trimmed temp FASTQs.
echo "[$TC_ID] done." >&2
