#!/bin/bash
# Build the combined bowtie2 index from <run_dir>/reference/combined_peptides.fasta
# (produced by build_reference.py; run_dir is configured in
# reference/paths.json, separate from this repo checkout). Pre-migration,
# no script in this repo did this step -- it was a manual/undocumented
# prerequisite; this closes that gap for the SLURM deployment.
#
# Submit with (PHIPSEQ_REPO_DIR is required -- see below):
#   export PHIPSEQ_REPO_DIR=/path/to/PhIPseq   # once per shell session
#   sbatch scripts/build_index.sbatch
# (after running `python scripts/build_reference.py` to produce the FASTA)
#
# PHIPSEQ_REPO_DIR is required -- SLURM copies a submitted script into a
# per-job spool directory on the compute node (e.g. /tmp/slurmd/...) and
# runs *that* copy, so the script's own location (${BASH_SOURCE[0]}) can't
# be used to find the repo checkout the way a plain bash script normally
# would.
#
# SLURM's own default output naming applies (this script has no static
# --output directive, since results/reference output now live in run_dir,
# not at a fixed path relative to the repo) -- pass --output= explicitly at
# submission time, or submit from within run_dir, if you want control over
# where the job's own stdout/stderr capture lands.
#
# #SBATCH lines are comments to bash but options to `sbatch` (see
# run_sample.sbatch for details).
# Name shown in `squeue`/`sacct` listings.
#SBATCH --job-name=phipseq_build_index
# CPU cores reserved; exported to the job as $SLURM_CPUS_PER_TASK.
#SBATCH --cpus-per-task=4
# Memory reserved for the job.
#SBATCH --mem=8G
# Wall-clock limit, hours:minutes:seconds (1 hour).
#SBATCH --time=01:00:00
# Strict mode: -e stop on the first failing command, -u error on unset
# variables, -o pipefail make a pipeline fail if any stage in it fails.
set -euo pipefail

# ${VAR:?message}: use PHIPSEQ_REPO_DIR, or stop with this error message if
# it is unset/empty.
PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"
# Read run_dir from reference/paths.json with a Python one-liner (-c "code")
# and capture its printed output with $(...). Uses whatever python3 is on
# PATH (stdlib json only, so no conda env needed). $PROJECT_DIR is expanded
# by bash into the Python string before Python runs.
RUN_DIR="$(python3 -c "import json; print(json.load(open('$PROJECT_DIR/reference/paths.json'))['run_dir'])")"

# Make bowtie2-build available, pinned to the version in
# reference/hpc_modules.json.
module load Bowtie2/2.5.4-GCC-14.2.0

# Create the index output folder (and parents) if it doesn't exist yet.
mkdir -p "$RUN_DIR/reference/bowtie2_index"
# Build the bowtie2 index (continued over 3 lines with trailing `\`):
#   --threads N    use the N CPUs SLURM allocated (--cpus-per-task above)
#   1st argument   the input FASTA of all peptide nucleotide sequences
#   2nd argument   output basename: writes combined.1.bt2 ... combined.rev.2.bt2
#                  (or .bt2l for very large references) into bowtie2_index/,
#                  which run_sample.sh then uses via `bowtie2 -x`.
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" \
  "$RUN_DIR/reference/combined_peptides.fasta" \
  "$RUN_DIR/reference/bowtie2_index/combined"
