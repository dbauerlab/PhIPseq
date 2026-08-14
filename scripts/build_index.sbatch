#!/bin/bash
# Build the combined bowtie2 index from reference/combined_peptides.fasta
# (produced by build_reference.py). Pre-migration, no script in this repo
# did this step -- it was a manual/undocumented prerequisite; this closes
# that gap for the SLURM deployment.
#
# Submit with: sbatch scripts/build_index.sbatch
# (after running `python scripts/build_reference.py` to produce the FASTA)
#SBATCH --job-name=phipseq_build_index
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=results/build_index_%j.slurm.log
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

module load Bowtie2/2.5.4-GCC-14.2.0

mkdir -p "$PROJECT_DIR/reference/bowtie2_index"
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" \
  "$PROJECT_DIR/reference/combined_peptides.fasta" \
  "$PROJECT_DIR/reference/bowtie2_index/combined"
