#!/bin/bash
# Differential enrichment of peptides vs the background library (see
# scripts/differential_enrichment.py and reference/de_params.json), per IP.
# Reads count_matrix.csv, so it must run after summarize.sbatch:
#
#   export PHIPSEQ_REPO_DIR=/path/to/PhIPseq   # once per shell session
#   sbatch --dependency=afterok:<summarize_job_id> scripts/differential.sbatch
#
# (or with no --dependency if count_matrix.csv is already up to date).
# PHIPSEQ_REPO_DIR is required for the same reason as in summarize.sbatch:
# SLURM runs a spooled copy of this script, so its own location can't be
# used to find the repo. More CPUs than summarize.sbatch: PyDESeq2 fits
# ~100k peptide-wise GLMs per IP and parallelises over SLURM_CPUS_PER_TASK.
#
# Output goes to <run_dir>/results/differential/ (per-contrast CSVs,
# differential_results.csv, enrichment_summary.csv, organism_summary.csv,
# figures and summary.md).

# #SBATCH lines are comments to bash but options to `sbatch` (see
# run_sample.sbatch for details).
# Name shown in `squeue`/`sacct` listings.
#SBATCH --job-name=phipseq_differential
# CPU cores reserved; exported as $SLURM_CPUS_PER_TASK, which
# differential_enrichment.py passes to PyDESeq2 as its worker count.
#SBATCH --cpus-per-task=8
# Memory reserved for the job.
#SBATCH --mem=32G
# Wall-clock limit, hours:minutes:seconds (4 hours).
#SBATCH --time=04:00:00

# Strict mode: -e stop on the first failing command, -u error on unset
# variables, -o pipefail make a pipeline fail if any stage in it fails.
set -euo pipefail

# ${VAR:?message}: use PHIPSEQ_REPO_DIR, or stop with this error message if
# it is unset/empty.
PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"

# Make the `conda` shell function available in this non-interactive job
# (`conda shell.bash hook` prints the setup code; eval runs it), then
# activate the phipseq env, which provides pydeseq2/pandas/matplotlib.
eval "$(conda shell.bash hook)"
conda activate phipseq

# Run the differential-enrichment analysis (settings from reference/de_params.json).
python "$PROJECT_DIR/scripts/differential_enrichment.py"
