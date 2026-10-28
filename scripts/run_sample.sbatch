#!/bin/bash
# SLURM array job: one task per row of reference/samples.tsv -- the sole
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
# Submit with:
#   export PHIPSEQ_REPO_DIR=/path/to/PhIPseq   # once per shell session
#   sbatch scripts/run_sample.sbatch
# (or: sbatch --export=ALL,PHIPSEQ_REPO_DIR=/path/to/PhIPseq scripts/run_sample.sbatch)
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
# capture lands. run_sample.sh's own cutadapt.log/bowtie2.log always land
# in run_dir/results/<tc_id>/ regardless.
#
# How #SBATCH lines work: to bash they are ordinary comments, but `sbatch`
# reads them (only while they appear before the first real command) as if
# they were command-line options. Options given on the sbatch command line
# override them, e.g. `sbatch --array=30-31 scripts/run_sample.sbatch`.
# The default job log is slurm-<jobid>_<taskid>.out in the submission dir.
#
# Name shown in `squeue`/`sacct` listings.
#SBATCH --job-name=phipseq_run_sample
# CPU cores reserved per array task (one task = one sample); exported to the
# job as $SLURM_CPUS_PER_TASK.
#SBATCH --cpus-per-task=4
# Memory reserved per task (the job is killed if it exceeds this).
#SBATCH --mem=32G
# Wall-clock limit per task, format days-hours:minutes:seconds (7 days).
#SBATCH --time=7-00:00:00
# Run one copy of this script per index in 1..29; each copy gets its own
# index in $SLURM_ARRAY_TASK_ID. NOTE: reference/samples.tsv currently has
# 31 sample rows, so with 1-29 the last two samples (PhIP_Seq1_35,
# PhIP_Seq1_36) are NOT run by a plain `sbatch` -- override with
# `--array=1-31` (or `--array=30-31`) at submission time.
#SBATCH --array=1-29
# Strict mode: -e stop on the first failing command, -u error on unset
# variables, -o pipefail make a pipeline fail if any stage in it fails.
set -euo pipefail

# ${VAR:?message} means "the value of VAR, but if it is unset or empty, print
# message to stderr and exit with an error" -- a required-variable check.
PROJECT_DIR="${PHIPSEQ_REPO_DIR:?Set PHIPSEQ_REPO_DIR to the repo checkout path before submitting (see header comment above)}"

# Pin exact tool versions from reference/hpc_modules.json -- keep this line
# in sync with that file by hand if the cluster's module versions change.
# `module load` (environment modules / Lmod) adjusts PATH and related
# environment variables so the named, pre-installed software versions
# (bowtie2, cutadapt, samtools) become runnable in this job. zcat/gzip/sort
# etc. used by run_sample.sh are standard system tools and need no module.
module load Bowtie2/2.5.4-GCC-14.2.0 cutadapt/4.2-GCCcore-11.3.0 SAMtools/1.22.1-GCC-14.2.0

# Activate the conda env built by scripts/setup_env.sh (provides the python3
# used for reference/*.json config loading in run_sample.sh). Adjust this if
# this cluster needs a different conda init step (e.g. `module load
# Anaconda3` before `conda activate` works).
# `conda shell.bash hook` prints the shell code that defines the `conda`
# shell function; `eval "$(...)"` runs that code here. This is needed in
# non-interactive batch jobs, where ~/.bashrc's conda init has not run.
eval "$(conda shell.bash hook)"
# Put the phipseq env's bin/ (python3 etc.) first on PATH.
conda activate phipseq

# SLURM_ARRAY_TASK_ID is 1-based; +1 skips samples.tsv's header line, so
# task 1 -> the first sample row, task 31 -> the last (samples.tsv currently
# has 31 sample rows; see the --array note above).
# $((...)) is bash integer arithmetic; LINE is the samples.tsv line number.
LINE=$((SLURM_ARRAY_TASK_ID + 1))
# Read that one line of samples.tsv into three variables:
#   sed -n "${LINE}p" FILE   print only line number LINE (-n = print nothing
#                            by default; `Np` = print line N)
#   < <(...)                 "process substitution": feed that command's
#                            output into `read` as if it were a file
#   read -r A B C            split the line on whitespace (tabs here) into the
#                            three columns tc_id, library, sample_prefix;
#                            -r = don't treat backslashes specially
read -r TC_ID LIBRARY SAMPLE_PREFIX < <(sed -n "${LINE}p" "$PROJECT_DIR/reference/samples.tsv")

# Log which sample this task is processing (goes to the SLURM job log).
echo "Task $SLURM_ARRAY_TASK_ID -> $TC_ID ($LIBRARY, $SAMPLE_PREFIX)"
# Run the per-sample pipeline (merge lanes -> cutadapt -> bowtie2 -> counts)
# on the full data (no 4th subsample argument).
bash "$PROJECT_DIR/scripts/run_sample.sh" "$TC_ID" "$LIBRARY" "$SAMPLE_PREFIX"
