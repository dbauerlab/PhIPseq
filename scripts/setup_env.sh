#!/bin/bash
# Build (or update) the "phipseq" conda environment from environment.yml.
# Run once on the cluster, and again whenever environment.yml changes.
#
# This only provisions the Python side of the pipeline (pandas, numpy,
# matplotlib, scipy, openpyxl, pydeseq2 -- see environment.yml).
# cutadapt/bowtie2/samtools are provided by
# the HPC module system instead -- see reference/hpc_modules.json and the
# `module load` lines in scripts/*.sbatch.
#
# Run it as: bash scripts/setup_env.sh  (needs the `conda` command on PATH).
#
# Strict mode: -e stop on the first failing command, -u error on unset
# variables, -o pipefail make a pipeline fail if any stage in it fails.
set -euo pipefail

# Absolute path of the repo checkout = parent of this script's scripts/
# folder (Python: Path(__file__).resolve().parent.parent). `$(...)` runs the
# inner command and substitutes its output; ${BASH_SOURCE[0]} is this
# script's path; dirname drops the filename; `cd .. && pwd` prints the
# parent's absolute path (inside a subshell, so our own cwd is unchanged).
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# `create` for a first-time build, `update` to bring an existing env in line
# with a changed environment.yml -- try create first, fall back to update.
#   -f FILE    the environment spec to build from
#   -n NAME    name of the env (overrides `name:` in the file; same value here)
#   a || b     run b only if a fails (a "fallback"; also stops `set -e` from
#              aborting the script when `create` fails because the env
#              already exists). The trailing `\` continues the command onto
#              the next line.
conda env create -f "$PROJECT_DIR/environment.yml" -n phipseq || \
  conda env update -f "$PROJECT_DIR/environment.yml" -n phipseq
