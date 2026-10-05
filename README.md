# Subsystems

Computational experiments for the paper *Per-State Signatures and Mean
Time-to-Failure: An Algorithm for Systems With Simultaneous Failures*
(G. Lagos, S. Ávila): signatures and mean time-to-failure (MTTF) at all 2^n
states of a monotone binary system, MTTF-threshold repair policies and their
long-run cost, under exchangeable failure models with a Markov survivor count
(Lévy-frailty Marshall-Olkin, LFMO).

The experiment design (systems, failure models and their calibration by the
pairwise beta factor, cost parameters) follows the companion paper on repair
policies with delays (repository `Simple-delayed`).

## Setup

    julia --project=. -e 'using Pkg; Pkg.instantiate()'

On the UAI HPC cluster, run the setup (and the tests) through `srun`, never on
the login node:

    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.test()'

Precompiling once before submitting an array job avoids every array task
compiling the package at the same time.

## Run tests

    julia --project=. -e 'using Pkg; Pkg.test()'

The tests compare the structure functions with an independent depth-first
search, the counts of working states with `test/reference/` (computed
independently of the Julia code, see `test/reference/W_reference.py`), the
minimal cuts with Table 9 of Yeh et al. (2024), and the calibration of the
failure models.

## Test systems

The 22 networks of Yeh et al. (2024) (`\cite{yeh2024new}`, Figure 7), as
listed in `data/systems/yeh2024_raw.txt`: two-terminal networks whose
**components are the edges**. Networks 6D and 7D are directed, the others
undirected; 6D is semi-coherent (4 irrelevant components); network 21 is
network 13 redrawn; network 22 is network 18 plus one edge. See
`data/systems/README.md`.

## Structural signatures of the complete systems

    sbatch templates/signatures.sbatch          # one array task per system, 1..22
    squeue -u $USER                             # wait until done
    srun --partition=compute-slim --ntasks=1 --mem=1G --time=00:10:00 \
        julia --project=. scripts/check_signatures.jl

Output: `data/signatures/<name>.csv` (committed, after copying them to a Mac),
with header comments and one row per `k = 0..n`: `W` (working states with `k`
failed components), `binomial`, `tail_signature` and `signature` (exact),
`signature_float`. `scripts/check_signatures.jl` compares every `W` with
`test/reference/yeh2024_W_reference.csv` and exits with code 1 on any
discrepancy.

## Failure models

`src/models.jl` (as in `Simple-delayed`): `Exp` and the LFMO models `LFMO-C`,
`LFMO-G`, `LFMO-IG`, `LFMO-S`, `BFR`, `LSh`, calibrated to Ψ(1) = 1 (unit mean
component lifetime) and pairwise beta factor β = 2 − Ψ(2); the experiments use
β ∈ {0.01, 0.05, 0.1, 0.2, 0.5}. Transition rates of the survivor-count chain
are computed in 512-bit BigFloat.

## Layout

- `src/` reusable code (systems, signatures, failure models)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/systems/` system definitions; `data/signatures/` computed signatures (committed)
- `test/reference/` independent reference values used by the tests
- `output/`, `logs/` results and job logs (git-ignored)
