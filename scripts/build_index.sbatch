#!/bin/bash
# Build the combined bowtie2 index from <run_dir>/reference/combined_peptides.fasta
# (produced by build_reference.py; run_dir is configured in
# reference/paths.json, separate from this repo checkout). Pre-migration,
# no script in this repo did this step -- it was a manual/undocumented
# prerequisite; this closes that gap for the SLURM deployment.
#
# Submit with: sbatch scripts/build_index.sbatch
# (after running `python scripts/build_reference.py` to produce the FASTA)
#
# SLURM's own default output naming applies (this script has no static
# --output directive, since results/reference output now live in run_dir,
# not at a fixed path relative to the repo) -- pass --output= explicitly at
# submission time, or submit from within run_dir, if you want control over
# where the job's own stdout/stderr capture lands.
#SBATCH --job-name=phipseq_build_index
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="$(python3 -c "import json; print(json.load(open('$PROJECT_DIR/reference/paths.json'))['run_dir'])")"

module load Bowtie2/2.5.4-GCC-14.2.0

mkdir -p "$RUN_DIR/reference/bowtie2_index"
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" \
  "$RUN_DIR/reference/combined_peptides.fasta" \
  "$RUN_DIR/reference/bowtie2_index/combined"
