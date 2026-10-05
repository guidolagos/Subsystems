#!/bin/bash
# Computational effort of Algorithm alg: one array task per system (scripts/benchmark_sweep.jl).
# Each task runs alone on one core: enumeration, sweep, and sweep with the MTTF of 31 models.
#
# Submit from the repository root:
#   sbatch --array=1-12,14-18,20,22 templates/benchmark_sweep.sbatch
#   sbatch --array=13,19,21 --mem=24G --time=12:00:00 templates/benchmark_sweep.sbatch
#SBATCH --job-name=bench-sweep
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
julia --project=. scripts/benchmark_sweep.jl "${SLURM_ARRAY_TASK_ID}"
echo "end: $(date -Is)"
