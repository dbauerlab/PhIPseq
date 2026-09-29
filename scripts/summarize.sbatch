#!/bin/bash
# Aggregate all samples' counts/logs into count_matrix.csv, qc_summary.csv,
# and figures. Light job: single-core, matplotlib over an already-aggregated
# matrix, no alignment/trimming -- resourced much smaller than
# run_sample.sbatch.
#
# Submit with (PHIPSEQ_REPO_DIR is required -- see below):
#   export PHIPSEQ_REPO_DIR=/path/to/PhIPseq   # once per shell session
#   sbatch --dependency=afterok:<run_sample_array_job_id> scripts/summarize.sbatch
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

#SBATCH --job-name=phipseq_summarize
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"

eval "$(conda shell.bash hook)"
conda activate phipseq

# These five don't depend on each other (summarize.py must run first to
# produce count_matrix.csv; the other four only read it), so running them
# sequentially in one job is simplest -- nothing in the pipeline suggests
# splitting them into separate parallel jobs would be worth the added
# complexity.
python "$PROJECT_DIR/scripts/summarize.py"
python "$PROJECT_DIR/scripts/completeness_figures.py"
python "$PROJECT_DIR/scripts/evenness_summary_figure.py"
python "$PROJECT_DIR/scripts/overview_figure.py"
python "$PROJECT_DIR/scripts/top_peptides.py"
