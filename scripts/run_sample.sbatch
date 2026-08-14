#!/bin/bash
# SLURM array job: one task per row of scripts/samples.tsv -- the sole
# orchestration mechanism for running all samples (there is no local/
# non-SLURM batch driver in this repo). SLURM's own array scheduling is
# left unthrottled (no `%N` on --array) -- the old MAX_PARALLEL=1
# sequential-only workaround existed only because the old SMB-mounted
# FASTQ share errored under concurrent reads, a limitation the
# HPC's storage doesn't share.
#
# Resources requested flat rather than tuned tightly: the 7-day wall-clock
# has negligible queue-time impact on this cluster, so there's no benefit
# to under-requesting. bowtie2 -p4 + cutadapt -j2 run sequentially within
# one sample (trim, then align), so 4 CPUs covers the peak.
#
# Submit with: sbatch scripts/run_sample.sbatch
#SBATCH --job-name=phipseq_run_sample
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --array=1-6
#SBATCH --output=results/run_sample_%A_%a.slurm.log
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Pin exact tool versions from reference/hpc_modules.json -- keep this line
# in sync with that file by hand if the cluster's module versions change.
module load Bowtie2/2.5.4-GCC-14.2.0 cutadapt/4.2-GCCcore-11.3.0 SAMtools/1.22.1-GCC-14.2.0

# Activate the conda env built by scripts/setup_env.sh (provides the python3
# used for reference/*.json config loading in run_sample.sh). Adjust this if
# this cluster needs a different conda init step (e.g. `module load
# Anaconda3` before `conda activate` works).
eval "$(conda shell.bash hook)"
conda activate phipseq

# SLURM_ARRAY_TASK_ID is 1-based; +1 skips samples.tsv's header line, so
# task 1 -> the first sample row, task 6 -> the last.
LINE=$((SLURM_ARRAY_TASK_ID + 1))
read -r TC_ID LIBRARY SAMPLE_PREFIX < <(sed -n "${LINE}p" "$PROJECT_DIR/scripts/samples.tsv")

echo "Task $SLURM_ARRAY_TASK_ID -> $TC_ID ($LIBRARY, $SAMPLE_PREFIX)"
bash "$PROJECT_DIR/scripts/run_sample.sh" "$TC_ID" "$LIBRARY" "$SAMPLE_PREFIX"
