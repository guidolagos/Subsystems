#!/bin/bash
# Structural signatures of the 22 systems of Yeh et al. (2024), by exhaustive
# enumeration of the 2^n states. One array task per system; task id = system id.
#
# Submit from the repository root:
#   sbatch templates/signatures.sbatch
# Re-submitting is safe: systems whose data/signatures/<name>.csv exists are skipped.
# To run only some systems:  sbatch --array=6,7 templates/signatures.sbatch
#
# Cost: 2^n connectivity checks. n <= 27 takes under half a minute; the three
# n = 30 systems (13, 19, 21) take 2-3 minutes each on one core.
#SBATCH --job-name=signatures
#SBATCH --partition=compute-slim
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out

module purge
module load softwares
module load Julia/1.10.4-linux-x86_64

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"

echo "host: $(hostname)  job: ${SLURM_JOB_ID}  task: ${SLURM_ARRAY_TASK_ID}  start: $(date -Is)"
julia --project=. --threads="${SLURM_CPUS_PER_TASK:-1}" \
    scripts/compute_signatures.jl "${SLURM_ARRAY_TASK_ID}"
echo "end: $(date -Is)"
