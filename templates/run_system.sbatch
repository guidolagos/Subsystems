#!/bin/bash
# MTTF-threshold policies: one array task per system (task id = system id, 1..22),
# each doing the 31 grid rows of its system (scripts/run_system.jl). The sweep of the
# states is cached in output/classes/<name>.bin; rows already done are skipped, so
# re-submitting is safe.
#
# Submit from the repository root, after precompiling once (see README):
#   sbatch --array=1-12,14-18,20,22 templates/run_system.sbatch              # n <= 27
#   sbatch --array=13,19,21 --mem=24G --time=12:00:00 templates/run_system.sbatch   # n = 30
# Memory: the two largest levels take ~9 GB for n = 30, ~1 GB for n = 27; the classes
# add to it (measure with sacct MaxRSS on a pilot).
#SBATCH --job-name=thresholds
#SBATCH --partition=compute-slim
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}  start: $(date -Is)"
julia --project=. scripts/run_system.jl "${SLURM_ARRAY_TASK_ID}" ${EXTRA:-}
echo "end: $(date -Is)"
