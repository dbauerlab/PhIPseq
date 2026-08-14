#!/bin/bash
# Build (or update) the "phipseq" conda environment from environment.yml.
# Run once on the cluster, and again whenever environment.yml changes.
#
# This only provisions the Python side of the pipeline (pandas, numpy,
# matplotlib, scipy, openpyxl). cutadapt/bowtie2/samtools are provided by
# the HPC module system instead -- see reference/hpc_modules.json and the
# `module load` lines in scripts/*.sbatch.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# `create` for a first-time build, `update` to bring an existing env in line
# with a changed environment.yml -- try create first, fall back to update.
conda env create -f "$PROJECT_DIR/environment.yml" -n phipseq || \
  conda env update -f "$PROJECT_DIR/environment.yml" -n phipseq
