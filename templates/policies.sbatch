#!/bin/bash
# Benchmark policies (scripts/run_policies.jl): rows dealt round-robin to NTASKS array
# tasks. Rows already done are skipped, so re-submitting is safe.
#
# Submit from the repository root, after precompiling once (see README):
#   main run, all 682 rows, state-based for n <= 23:
#     sbatch templates/policies.sbatch
#   pilot (task 1 only, a few rows):
#     sbatch --array=1 --export=ALL,EXTRA="--rows=1:40" templates/policies.sbatch
#   state-based for n = 24, 26, 27 (systems 20, 14, 17), one row per task, after a pilot:
#     sbatch --array=1-93%30 --mem=16G --time=12:00:00 \
#       --export=ALL,NTASKS=93,EXTRA="--min-n=24 --max-n=27 --state-max-n=27 --out=large" templates/policies.sbatch
# The --array range must be 1-NTASKS.
#SBATCH --job-name=policies
#SBATCH --partition=compute-slim
#SBATCH --array=1-40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out

NTASKS=${NTASKS:-40}

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}/${NTASKS}  start: $(date -Is)"
julia --project=. scripts/run_policies.jl "${SLURM_ARRAY_TASK_ID}" "${NTASKS}" ${EXTRA:-}
echo "end: $(date -Is)"
