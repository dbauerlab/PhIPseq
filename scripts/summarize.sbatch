#!/bin/bash
# Aggregate all samples' counts/logs into count_matrix.csv, qc_summary.csv,
# and figures. Light job: single-core, matplotlib over an already-aggregated
# matrix, no alignment/trimming -- resourced much smaller than
# run_sample.sbatch.
#
# Submit with (PHIPSEQ_REPO_DIR is required -- see below):
#   export PHIPSEQ_REPO_DIR=/path/to/PhIPseq   # once per shell session
#   sbatch --dependency=afterok:<run_sample_array_job_id> scripts/summarize.sbatch
# (--dependency=afterok:<id> makes SLURM hold this job until job <id> -- for
# an array job, every one of its tasks -- has finished successfully.)
#
# PHIPSEQ_REPO_DIR is required -- SLURM copies a submitted script into a
# per-job spool directory on the compute node (e.g. /tmp/slurmd/...) and
# runs *that* copy, so the script's own location (${BASH_SOURCE[0]}) can't
# be used to find the repo checkout the way a plain bash script normally
# would.
#
# SLURM's own default output naming applies (no static --output directive
# here, since results/ now lives in run_dir -- configured in
# reference/paths.json -- not at a fixed path relative to the repo) --
# pass --output= explicitly at submission time, or submit from within
# run_dir, if you want control over where the job's own stdout/stderr
# capture lands.

# #SBATCH lines are comments to bash but options to `sbatch` (see
# run_sample.sbatch for details).
# Name shown in `squeue`/`sacct` listings.
#SBATCH --job-name=phipseq_summarize
# A single CPU core (the scripts below are single-threaded).
#SBATCH --cpus-per-task=1
# Memory reserved for the job.
#SBATCH --mem=8G
# Wall-clock limit, hours:minutes:seconds (2 hours).
#SBATCH --time=02:00:00

# Strict mode: -e stop on the first failing command (so a failing script
# stops the ones after it), -u error on unset variables, -o pipefail make a
# pipeline fail if any stage in it fails.
set -euo pipefail

# ${VAR:?message}: use PHIPSEQ_REPO_DIR, or stop with this error message if
# it is unset/empty.
PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"

# Make the `conda` shell function available in this non-interactive job
# (`conda shell.bash hook` prints the setup code; eval runs it), then
# activate the phipseq env so `python` below has pandas/matplotlib/etc.
eval "$(conda shell.bash hook)"
conda activate phipseq

# These five run in a fixed order: summarize.py must run first to produce
# count_matrix.csv (and count_matrix_filtered.csv); the other four only read
# its outputs and don't depend on each other, so running them
# sequentially in one job is simplest -- nothing in the pipeline suggests
# splitting them into separate parallel jobs would be worth the added
# complexity.
# Builds count_matrix.csv, count_matrix_filtered.csv, filter_report.csv,
# qc_summary.csv, summary.md and the core figures (incl. PCA).
python "$PROJECT_DIR/scripts/summarize.py"
# Coverage / rarefaction / detection-histogram figures.
python "$PROJECT_DIR/scripts/completeness_figures.py"
# Per-sample evenness figure.
python "$PROJECT_DIR/scripts/evenness_summary_figure.py"
# Combined coverage + evenness overview figure.
python "$PROJECT_DIR/scripts/overview_figure.py"
# Top-variable-peptide heatmaps + top-abundant-peptide tables
# (top_abundant_peptides.csv, top_abundant_peptides_wide.csv).
python "$PROJECT_DIR/scripts/top_peptides.py"
