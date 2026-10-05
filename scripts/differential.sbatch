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

#SBATCH --job-name=phipseq_differential
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"

eval "$(conda shell.bash hook)"
conda activate phipseq

python "$PROJECT_DIR/scripts/differential_enrichment.py"
