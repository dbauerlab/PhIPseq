# SLURM/HPC Migration Plan

Working notes for the `nemo` branch: review the current macOS/local pipeline
(documented in [README.md](README.md)) and redeploy it on a SLURM-scheduled
HPC cluster. This is a plan only — no code has been changed yet.

## Goals

1. Reproducible conda environment instead of an untracked local `.venv`.
2. Retire scripts that only existed to work around macOS/local-storage
   quirks or one-off incidents, now that adapter detection is being replaced
   by a user-maintained `reference/adapters.json`.
3. Move all hardcoded paths/thresholds/constants out of script bodies and
   into JSON config under `reference/`.
4. Convert the bash driver scripts to SLURM submission scripts.
5. Merge sequencing lanes before trimming/alignment (rather than
   processing each lane separately and combining counts afterward), and
   align read pairs jointly in paired-end mode (rather than trimming and
   aligning R1/R2 independently as single-end reads).
6. Comment every migrated/rewritten script to explain what the code is
   doing, throughout.

## 1. Conda environment

**Add `environment.yml`** at the repo root for the Python side only —
`cutadapt`/`bowtie2`/`samtools` are provided by the HPC module system
instead (see §1a), not conda, since specific pinned module builds are
already available on the cluster:

```yaml
name: phipseq
channels:
  - conda-forge
dependencies:
  - python              # unpinned: always resolve the current stable conda-forge release
  - pandas>=2.2
  - numpy>=2.1
  - matplotlib>=3.9
  - scipy>=1.14
  - openpyxl>=3.1        # pandas.read_excel engine for build_reference.py
```

- Python version deliberately left unpinned, per instruction to always take
  the most up-to-date stable release available at env-build time, rather
  than freezing a specific minor version now.
- Package floors above are recommended current-stable versions (as of this
  plan being written) compatible with an unpinned recent Python 3 — not
  confirmed against what the existing `results/` were generated with (no
  version pins exist anywhere in the current repo). Treat as a reasonable
  starting point; bump the floors if `conda`'s solver pulls in anything
  older, and re-check before final pinning if exact reproducibility of the
  already-generated `results/` matters.
- `bioconda`/`cutadapt`/`bowtie2`/`samtools` removed from the env — these
  three come from `module load` on the HPC instead (§1a), so the conda env
  only needs to satisfy the pure-Python scripts.

**1a. HPC module versions — `reference/hpc_modules.json`**

Pin the exact module versions confirmed available on the cluster's module
system, as the single source of truth for every `.sbatch` script's
`module load` lines:

```json
{
  "bowtie2": "Bowtie2/2.5.4-GCC-14.2.0",
  "cutadapt": "cutadapt/4.2-GCCcore-11.3.0",
  "samtools": "SAMtools/1.22.1-GCC-14.2.0"
}
```

Each `.sbatch` script loads these three via `module load`, e.g.:

```bash
module load Bowtie2/2.5.4-GCC-14.2.0 cutadapt/4.2-GCCcore-11.3.0 SAMtools/1.22.1-GCC-14.2.0
```

(hardcoded in the `.sbatch` header to match `reference/hpc_modules.json` —
if the cluster's module versions change, update both together; a helper
that reads the JSON and emits `module load` lines via `jq` is an option if
keeping the two in sync by hand becomes error-prone).

**Add `scripts/setup_env.sh`** — thin wrapper to build the env
non-interactively on a login/build node:

```bash
#!/bin/bash
set -euo pipefail
conda env create -f environment.yml -n phipseq || \
  conda env update -f environment.yml -n phipseq
```

- Replaces the current `PROJECT_DIR/.venv` convention referenced directly in
  `run_sample.sh`. Every script that currently shells out to
  `"$PROJECT_DIR/.venv/bin/python"` needs that hardcoded path replaced with
  `conda run -n phipseq python` (or an activated env on the compute node —
  decide per SLURM script, see §4).

## 2. Scripts to retire

| Script | Disposition | Reason |
|---|---|---|
| `scripts/detect_adapters.py` | **Delete** | Adapter sequences will be entered directly into `reference/adapters.json` by hand; no script needs to (re)derive them. Confirmed: delete outright rather than keep unused. |
| `scripts/wait_and_run.sh` | **Remove** | One-off recovery wrapper for a specific historical incident (an orphaned `TC_1` run from an earlier driver bug — see its own header comment). Not a general entry point; `run_all.sh`/its SLURM replacement is sufficient for a normal run. |

Both are currently documented as script entries in [README.md](README.md#detect_adapterspy)
— once removed, that documentation (and the Mermaid data-flow diagram) needs
a follow-up README update.

Also resolved: `run_sample.sh`'s `cat | gzcat` two-stage pipe exists
specifically to dodge a macOS/SMB `close()`/EBADF bug (per its comment).
The HPC has `zcat` available, so this migration replaces the workaround
with plain `zcat` (or `zcat "$fq"` directly in `read_fastq()`) — the
`cat |` indirection and its associated comment can be dropped along with
it.

## 3. Config consolidation

Everything currently hardcoded as a module-level Python constant or shell
variable moves into JSON under `reference/`, one settings file per concern:

- **`reference/paths.json`** — external input locations currently hardcoded
  in-script (FASTQ directory in `run_sample.sh`; the CoV `.xlsx` and Vir3
  `.csv` source reference files in `build_reference.py`). Not reproducing
  the literal current macOS/Dropbox paths here since they're lab-storage-
  specific and will change on the new cluster anyway. Confirmed: ship this
  file with placeholder values (e.g. `"fastq_dir": "TODO"`) for now — the
  user will fill in the real HPC paths at the point they actually run the
  pipeline, not during this migration.
- **`reference/run_params.json`** — sequencing/alignment tunables currently
  hardcoded in `run_sample.sh` / `run_all.sh`: `LANES` (`L007`, `L008`),
  `BOWTIE_THREADS`, `CUTADAPT_THREADS`. `MAX_PARALLEL` and the per-lane
  retry count/backoff (currently 5 attempts / 20s) are both dropped
  entirely rather than migrated — SLURM's own array scheduler and job-level
  failure handling take over both roles (§4, §5a, §8).
- **`reference/build_reference_params.json`** — the CoV/Vir3 vector anchor
  sequences (`COV_ANCHOR_5`/`COV_ANCHOR_3`/`VIR3_ANCHOR_5`/`VIR3_ANCHOR_3`
  in `build_reference.py`) and the minimum Vir3 insert length (currently 15nt,
  inline in the code) — these are library-design constants, not per-run
  parameters, but still shouldn't live inside the script body.
- **`reference/figure_params.json`** — figure thresholds/styling currently
  duplicated across `evenness_summary_figure.py` and `overview_figure.py`:
  `OUTLIER_THRESHOLD` (currently 10, defined separately in each of the two
  scripts — worth deduplicating into one config value while migrating) and
  `LIBRARY_COLORS` (`{"CoV": "tab:blue", "Vir3": "tab:orange"}`).
- **`reference/adapters.json`** — already exists, already JSON, no change in
  location; just changes from generated (by `detect_adapters.py`) to
  hand-maintained per the user's stated plan.

Each script currently reading these as constants needs a small loader
(`json.load(open(PROJECT_DIR / "reference" / "....json"))`) added at the top
in place of the literal constant — mechanical, one script at a time.

## 4. SLURM conversion

Current shape: `run_all.sh` loops over `scripts/samples.tsv` and launches
`run_sample.sh` in the background, one sample at a time (`MAX_PARALLEL=1`,
forced by the old SMB mount's concurrency limits, not a compute constraint).
Proposed SLURM shape:

- **`scripts/run_sample.sbatch`** — `run_sample.sh`'s per-sample body,
  wrapped with an `#SBATCH` header. Resources per user decision (queue-time
  impact of the 7-day wall-clock is negligible on this cluster, so no need
  to under-request):
  ```bash
  #SBATCH --job-name=phipseq_%a
  #SBATCH --cpus-per-task=4
  #SBATCH --mem=32G
  #SBATCH --time=7-00:00:00
  #SBATCH --array=1-6         # one task per samples.tsv row; no throttle (see below)
  ```
  Sample row lookup by `SLURM_ARRAY_TASK_ID` replaces `run_all.sh`'s bash
  loop. `MAX_PARALLEL` is not carried forward at all — `--array=1-6` is left
  unthrottled (no `%N`) and the scheduler is free to run all 6 tasks
  concurrently if resources allow.
- **`scripts/summarize.sbatch`** — `summarize.py` +
  `completeness_figures.py` + `evenness_summary_figure.py` +
  `overview_figure.py`, submitted with a job dependency on the array job:
  ```bash
  sbatch --dependency=afterok:<run_sample_array_job_id> scripts/summarize.sbatch
  ```
  These four scripts don't depend on each other (per
  [README.md § Script running order](README.md#script-running-order)) and
  are all short/single-core (matplotlib figure generation over an
  already-aggregated matrix), so one job running them sequentially is
  simplest; splitting into 4 parallel job steps is possible later if this
  becomes a bottleneck, but nothing in the current code suggests it would.
- **Bowtie2 index build** — currently unresolved/manual (no script in the
  repo calls `bowtie2-build`, confirmed in README). This migration is the
  right point to add an explicit `scripts/build_index.sbatch` (or fold into
  `build_reference.py`'s job) rather than carrying the gap forward.
- Replace the `.venv/bin/python` calls with the conda env from §1, e.g.
  `conda run -n phipseq python scripts/summarize.py` or activating the env
  at the top of each `.sbatch` script.

## 5. Read pre-processing changes: lane merging + paired-end alignment

Current `run_sample.sh` logic (per mate, R1 and R2 handled independently):
for each lane (`L007`, `L008`) separately — decompress → `cutadapt -g
<adapter> --discard-untrimmed` → `bowtie2 -U -` (single-end, `--norc` for
R1 / `--nofw` for R2 to fake strand-specificity) → `samtools view -F 4` →
`cut -f3`; then concatenate the two lanes' hit lists and `sort | uniq -c`
into `R1_counts.tsv` / `R2_counts.tsv`. R1 and R2 are two entirely separate
alignment runs against the same index, never treated as mate pairs.

Two changes requested, both in `run_sample.sh`:

**5a. Merge lanes before processing, not after.** Concatenate `L007` +
`L008` per mate *before* trimming/alignment (gzip-concatenation is valid —
`zcat L007.fastq.gz L008.fastq.gz` decompresses both as one stream, or the
two `.gz` files can be `cat`-concatenated directly and decompressed once),
rather than running the full trim→align→count pipeline once per lane and
combining count files afterward. This replaces the current per-lane `for
lane in "${LANES[@]}"` loop + hit-list concatenation with a single merged
stream per mate. The per-lane retry loop (5 attempts / 20s backoff) is
dropped entirely rather than migrated to wrap the merged stream — see §8.

Temp storage is still needed, though — see §5c. `reference/run_params.json`
(§3) still records `LANES`, but it now drives which raw FASTQ files get
merged, not a per-lane retry loop.

**5b. Align in paired-end mode.** Replace the two independent single-end
alignments (R1 with `--norc`, R2 with `--nofw`) with one paired-end
`bowtie2` invocation using `-1 <trimmed R1> -2 <trimmed R2>`, so mate pairs
are aligned jointly rather than each mate scored independently against the
index. Confirmed/resolved:
- Trimming (`cutadapt`) runs in paired mode (`cutadapt -g <R1 adapter> -G
  <R2 adapter> --discard-untrimmed -o <tmp R1> -p <tmp R2> <merged R1>
  <merged R2>`), trimming both mates together in one invocation instead of
  today's two separate single-mate calls.
- Paired-mode `cutadapt`'s default pair-filtering behavior (discard a pair
  if *either* mate fails a filter, e.g. `--discard-untrimmed`) is accepted
  as-is — no need to override with `--pair-filter=both`.
- The `--norc`/`--nofw` single-end orientation hack is dropped — confirmed
  the library prep is standard mate1-forward/mate2-reverse, i.e. `bowtie2`
  paired mode's default `--fr` orientation is correct as-is, no override
  needed.
- **Downstream impact on count matrices, still to be implemented, tracked
  but out of scope for this document**: today's output is two separate
  files/matrices, `R1_counts.tsv`/`R2_counts.tsv` →
  `count_matrix_R1.csv`/`count_matrix_R2.csv`, and `summarize.py`'s
  downstream figure scripts all key off the R1 matrix specifically (per
  [README.md § Shared inputs/outputs](README.md#shared-inputsoutputs)).
  Paired-end alignment produces one alignment (and so one count) per read
  pair/fragment, not one per mate — so `run_sample.sh`'s output becomes a
  single `results/<tc_id>/counts.tsv` → `results/count_matrix.csv` (see
  §5c for the counting-logic change this requires), which ripples into
  `completeness_figures.py`, `evenness_summary_figure.py`, and
  `overview_figure.py`, all of which currently assume an R1/R2 split.
  That redesign needs its own pass through those three figure scripts when
  implemented — the parsing side of the same ripple, into `summarize.py`
  specifically, is planned in §5d.

**5c. Temp files for the `cutadapt` → `bowtie2` handoff, and temp-directory
approach.** Today's per-lane pipeline is a single shell pipe
(`cutadapt | bowtie2 | samtools`) because single-end `bowtie2 -U -` reads
one stdin stream. Paired-mode `bowtie2 -1 <file1> -2 <file2>` needs two
separate, complete inputs, so a single pipe no longer works. Resolved:
write `cutadapt`'s paired output to two temp FASTQ files (`-o`/`-p`, not
stdout), then invoke `bowtie2 -1 <tmp R1> -2 <tmp R2>` reading from those
files — simpler and easier to debug than named pipes/process substitution,
at the cost of temp FASTQ files briefly touching disk instead of a pure
in-memory pipe.
- This still needs a temp directory (the old `TMP_DIR=/tmp/phipseq_tmp/$TC_ID`
  isn't going away, just changing shape): drop the hardcoded fixed path and
  create a fresh, unique directory per job with `mktemp -d` instead — it
  automatically honors `$TMPDIR` if the scheduler sets one, falls back to
  `/tmp` otherwise, avoids collisions between array tasks/reruns without
  guessing at cluster-specific scratch conventions, and should be removed
  with a `trap ... EXIT` cleanup rather than a manual `rm -f` at the end of
  the function (so it's cleaned up even if the script exits early on a
  failure).

**5d. Counting logic — one count per read pair, not per alignment record.**
`samtools view -F 4 | cut -f3` on a paired-mode BAM emits one row per
*mate* alignment (both R1 and R2 records for every aligned pair), so the
current `sort | uniq -c` would double-count every fragment relative to
before. Resolved approach: filter to one record per pair before counting —
`samtools view -F 4 -f 64` (`-f 64`/`0x40` = first-in-pair) selects exactly
one row per aligned pair, keeping the existing `cut -f3 | sort | uniq -c`
shape downstream unchanged. (`-F 4` still excludes pairs where the R1 mate
itself didn't align; whether to also require the *pair* to align
concordantly, e.g. via `-f 2`, is a implementation-time detail to confirm
against what "detected" should mean for the QC report — flagged for
attention during implementation, not a blocker to this plan.)

**5e. `summarize.py` log/count-matrix parsing update.**
Confirmed as in-scope for this plan (design only — implementation happens
alongside §5b–§5d). Current `summarize.py` structure (verified in the
script): `sample_readstats(tc_id, mate)` loops `for lane in LANES`, parsing
`cutadapt_{mate}_{lane}.log` and `bowtie2_{mate}_{lane}.log` per lane and
summing across the two lanes, called once for `"R1"` and once for `"R2"`;
`qc_summary.csv` then carries `raw_reads_R1`, `trim_rate_R1`,
`map_rate_R1`, `multimap_rate_R1`, and `map_rate_R2` as separate columns,
and all peptide-count QC (`n_peptides_detected`, `pct_library_covered`,
`n_dropouts`, `median_count_detected`, `gini_own_library`,
`pct_mapped_to_other_library`) is computed from `mat_r1` only.

Post-migration, per sample there is one merged, paired-mode `cutadapt.log`
and one paired-mode `bowtie2.log` (no per-lane, no per-mate split — see
§5a/§5b), and one `counts.tsv`/`count_matrix.csv` (no R1/R2 split — see
§5d). This requires:
- Deleting the `for lane in LANES` summation loop entirely — one log file
  per sample per tool, nothing to sum.
- Replacing `parse_cutadapt_log` with a paired-mode parser: cutadapt's
  paired-mode summary reports pair-level fields (`Total read pairs
  processed`, `Pairs written (passing filters)`) rather than the current
  single-end field names (`Total reads processed`, `Reads written`) —
  different regex targets, same purpose (raw vs. trimmed count).
- Replacing `parse_bowtie2_log` with a paired-mode parser: paired `bowtie2`
  stdout reports concordant/discordant pair-alignment rates (and a
  secondary per-mate breakdown for pairs that didn't align concordantly),
  structurally different from the current single-end `aligned 0
  times`/`exactly 1 time`/`>1 times` breakdown — the exact fields to carry
  into `qc_summary.csv` (e.g. a `concordant_map_rate` replacing today's
  `map_rate_R1`/`map_rate_R2` pair) need to be picked during implementation
  from a real paired-mode `bowtie2` log, not guessed here.
- `qc_summary.csv`'s schema changes from the current R1/R2-suffixed columns
  to one set of pair-level columns (naming TBD at implementation time,
  e.g. `raw_read_pairs`/`trim_rate`/`concordant_map_rate` in place of
  `raw_reads_R1`/`trim_rate_R1`/`map_rate_R1`/`multimap_rate_R1`/
  `map_rate_R2`) — this is a schema change `results/summary.md`'s QC table
  needs to follow.

## 6. Code documentation standard

Every script touched during this migration (new `.sbatch` files, rewritten
`run_sample.sh`, any script edited for config-loading per §3, or built-out
data-flow per §5) must be commented throughout to explain what the code is
doing, not left comment-sparse as several of the current scripts are.

## 7. Suggested execution order (post-migration)

1. `scripts/setup_env.sh` (once, or whenever `environment.yml` changes)
2. `python scripts/build_reference.py` → `reference/combined_peptides.fasta`,
   `reference/combined_metadata.csv`
3. `sbatch scripts/build_index.sbatch` (new) → `reference/bowtie2_index/combined`
4. Hand-maintained `reference/adapters.json` (no script — confirm values are
   present before step 5)
5. `sbatch scripts/run_sample.sbatch` (array job over `scripts/samples.tsv`)
6. `sbatch --dependency=afterok:<job_id> scripts/summarize.sbatch`

## 8. Resolved decisions

- **Concurrency**: HPC storage doesn't have the old SMB share's
  concurrent-read limitation, so `MAX_PARALLEL` is dropped and
  `--array=1-6` runs unthrottled.
- **Decompression**: `zcat` is available on the HPC; the `cat | gzcat`
  macOS/SMB workaround is dropped in favor of plain `zcat`.
- **Resource requests**: `--cpus-per-task=4`, `--mem=32G`,
  `--time=7-00:00:00` for `run_sample.sbatch` (7-day wall-clock has
  negligible queue-time impact on this cluster, so it's requested flat
  rather than tuned tightly).
- **`detect_adapters.py`**: delete outright (§2).
- **`cutadapt`/`bowtie2`/`samtools` versions**: pinned via HPC modules, not
  conda — see `reference/hpc_modules.json` (§1a).
- **Python package versions**: no version pins exist anywhere in the
  current repo to reconcile against, so the recommended floors from §1
  stand as final: `python` unpinned (always take the current stable
  release); `pandas>=2.2`, `numpy>=2.1`, `matplotlib>=3.9`, `scipy>=1.14`,
  `openpyxl>=3.1`.
- **`summarize.sbatch` resources**: confirmed — the smaller placeholder
  applies, since `summarize.py` + the three figure scripts are light
  (single-core, matplotlib over an already-aggregated matrix, no
  alignment/trimming):
  ```bash
  #SBATCH --job-name=phipseq_summarize
  #SBATCH --cpus-per-task=1
  #SBATCH --mem=8G
  #SBATCH --time=02:00:00
  ```
- **Per-lane retry loop removed**: agreed with the user's assessment — the
  5-attempts/20s-backoff retry in `run_sample.sh` exists specifically to
  paper over a transient SMB-mount `close()`/EBADF bug on the old
  Dropbox/macOS setup (per its own comment), not a general-purpose
  alignment safety net. That failure mode is tied to the old network
  share, not expected on the HPC's storage, and SLURM already provides the
  idiomatic replacement for job-level failures: a failed array task shows
  a nonzero exit in `sacct`/`seff`, and can be resubmitted individually
  (or via `--requeue`) rather than looping inside the script. Dropped
  entirely rather than migrated (§5a).
- **Real path values**: `reference/paths.json` (and any other config
  needing the actual FASTQ/source-reference locations) ships with
  placeholders; the user fills in real HPC paths later, at actual run time
  (§3).
- **Temp directory**: still needed — paired-mode `cutadapt` output now
  lands in temp FASTQ files that `bowtie2 -1/-2` reads from (§5c) — but
  the fixed `/tmp/phipseq_tmp/$TC_ID` path is replaced with a per-job
  `mktemp -d` (honors `$TMPDIR` if the scheduler sets one, else `/tmp`),
  cleaned up via `trap ... EXIT` (§5c).
- **`cutadapt`/`bowtie2` handoff**: temp files (`cutadapt -o/-p` to two
  temp FASTQs, then `bowtie2 -1/-2` reads them), not named pipes/process
  substitution — simpler to implement and debug (§5c).
- **Counting logic**: adjusted to select one BAM record per aligned pair
  (`samtools view -F 4 -f 64`, first-in-pair) before the existing
  `cut -f3 | sort | uniq -c`, to avoid double-counting each fragment now
  that both mates' alignment records appear in the same BAM (§5d).
- **`cutadapt` paired-mode default pair-filtering** (`--pair-filter=any` —
  discard a pair if either mate fails a filter): accepted as-is, no
  override needed (§5b).
- **`bowtie2` paired-mode default orientation** (`--fr`, mate1
  forward/mate2 reverse): confirmed correct for this library prep, no
  override needed (§5b).
- **`summarize.py` parsing update**: in scope for this plan at the design
  level — drop the per-lane/per-mate summation entirely, write new
  paired-mode `cutadapt`/`bowtie2` log parsers, and change `qc_summary.csv`
  to one set of pair-level columns instead of the current R1/R2-suffixed
  columns (§5e).

## 9. Repo checkout vs. run directory

Implemented post-migration, at the user's request: the repo checkout and
the pipeline's run output are kept in separate directories, rather than
everything living under wherever the repo is cloned.

- **Mechanism**: a new `run_dir` key in `reference/paths.json` (placeholder,
  same as the other path values in that file — see §3/§8). Every Python
  script and `run_sample.sh` resolves `run_dir` from this file (found via
  the repo checkout's own location, `Path(__file__).resolve().parent.parent`),
  then writes/reads all generated data there instead of relative to the
  script's own location.
- **What moved out of the repo**: `results/` (per-sample logs/counts,
  `count_matrix.csv`, `qc_summary.csv`, `summary.md`, figures) plus the
  *generated* reference data — `combined_peptides.fasta`, the bowtie2
  index, and `combined_metadata.csv` (previously tracked in git despite
  being generated — that convention ends here). Hand-maintained config
  (`reference/*.json` other than the generated files, `adapters.json`,
  `scripts/samples.tsv`) stays in the repo checkout, since it's small,
  version-controlled config rather than run output.
- **`.sbatch` `--output` directives dropped**: `run_sample.sbatch`,
  `summarize.sbatch`, and `build_index.sbatch` no longer set a static
  `#SBATCH --output=results/...` path, since `results/` is no longer at a
  fixed location relative to the repo (it's wherever `run_dir` points, only
  known once the script parses JSON at runtime — too late for an `#SBATCH`
  directive, which SLURM parses before any script code runs). SLURM's own
  default output naming applies instead; pass `--output=` explicitly at
  submission time, or submit from within `run_dir`, for control over where
  the job log lands.
- **`build_index.sbatch`** needed its own small `python3 -c` call to resolve
  `run_dir`, since (unlike `run_sample.sbatch`/`summarize.sbatch`) it builds
  paths itself rather than delegating entirely to another script that
  already does the JSON loading.
- **Pre-migration `results/`/`reference/` content already checked into this
  repo is untouched** by this change — it predates the run_dir split
  entirely and is called out as a legacy/historical artifact in
  [README.md](README.md#overview), not deleted or moved.
