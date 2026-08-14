#!/bin/bash
# Aggregate all samples' counts/logs into count_matrix.csv, qc_summary.csv,
# and figures. Light job: single-core, matplotlib over an already-aggregated
# matrix, no alignment/trimming -- resourced much smaller than
# run_sample.sbatch.
#
# Submit with a dependency on the run_sample.sbatch array job, so it only
# starts once every sample has finished:
#   sbatch --dependency=afterok:<run_sample_array_job_id> scripts/summarize.sbatch
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

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

eval "$(conda shell.bash hook)"
conda activate phipseq

# These four don't depend on each other (summarize.py must run first to
# produce count_matrix.csv; the other three only read it), so running them
# sequentially in one job is simplest -- nothing in the pipeline suggests
# splitting them into separate parallel jobs would be worth the added
# complexity.
python "$PROJECT_DIR/scripts/summarize.py"
python "$PROJECT_DIR/scripts/completeness_figures.py"
python "$PROJECT_DIR/scripts/evenness_summary_figure.py"
python "$PROJECT_DIR/scripts/overview_figure.py"
