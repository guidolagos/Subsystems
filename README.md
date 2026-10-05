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

## Per-state signatures and MTTF (Algorithm alg)

`src/sweep.jl` sweeps the 2^n states by number of working components and computes
the counting vector W^x of every state from those of its sons (Lemma lemma:Wk),
evaluating Φ once per state and keeping two levels in memory (UInt16/UInt32
matrices, one column per state; ~9 GB for n = 30). Working states are grouped
into classes of equal counting vector, with multiplicities: under an LFMO model
the MTTF of a state, and its contribution to the cost of the MTTF-threshold
policies, depend on the state only through its class. A state defines a minimal
path set iff W^x_0 = 1 and W^x_1 = 0. Classes are cached in
`output/classes/<name>.bin`; `output/sweep/<name>.csv` has the run time, the number
of classes and of minimal path sets.

## MTTF-threshold policies

`src/thresholds.jl`: for an LFMO model, m(x) = Σ_k W^x_k/binomial(w,k) δ^w_k with
δ^w_k = E T_{k+1:w} − E T_{k:w} computed from the survivor-count chain (nonnegative
terms). Since m is monotone, the continuation set {m(x) ≥ γ m(1)} of π^γ is an
up-set, and each of its states is occupied as in the chain without repairs. This
gives, in one pass over the classes sorted by m, the mean cycle length T, the mean
number of failed components at the repair N and the probability P that the repair
is corrective, at every breakpoint γ of γ ↦ LTMC^γ (the distinct values m(x)/m(1)),
and

    LTMC^γ = (c_set + c_cmp N + c_sys P) / T,   c_sys = 1, c_cmp = ratio/n, c_set = setup,

for the cost pairs ratio ∈ `COST_RATIOS` (1/10, 1/5, 1/2, 2/3, 1, 4/3, 2, 4, 10)
and setup ∈ `COST_SETUPS` (0, 0.05, 0.1, 0.2, 0.5), as in the companion paper.

The grid `data/grid.csv` (committed; `scripts/make_grid.jl`) has one row per
system (22) and model (31: `Exp`, and `LFMO-C`, `LFMO-G`, `LFMO-IG`, `LFMO-S`, `BFR`,
`LSh` with β ∈ {0.01, 0.05, 0.1, 0.2, 0.5}): 682 rows; system k owns rows
31(k−1)+1..31k. Each row writes to `output/thresholds/`:

- `row_<id>.csv`: one line per cost pair: `n_breakpoints`, the best breakpoint
  (`best_index`, `best_gamma`; any γ in (`best_gamma_lo`, `best_gamma`] gives the same
  policy), `ltmc_best`, `ltmc_gamma1` (γ = 1: repair at the first relevant failure),
  `ltmc_corrective` (γ → 0+), and `local_minima` of LTMC over the breakpoints (1 = unimodal);
- `row_<id>_curve.csv`: T, N, P at γ = 0, 0.005, ..., 1 (γ = 0 means 0+), from which
  the curve of any cost pair can be drawn;
- `row_<id>_breaks.csv`: all breakpoints, when there are at most 10,000.

Files carry `# thresholds_version`; readers refuse other versions.

Validation: `test/reference/thresholds_reference.py` follows the paper literally
(MTTF by the recursion eq:mttf, without signatures; T, N, P from the linear systems
sys:m rep and sys:c rep), with rates in 60-digit arithmetic; `Pkg.test()` compares
24 cases (systems 1, 3, 6D, 9; all seven models) with relative tolerance 1e-10.
`test/reference/thresholds_port.py` is a Float64 Python port of the Julia design
(classes, occupation times), checked against the reference to 5e-13. The engine
also checks, for every row, that T(0+) = m(1), P(0+) = 1 and that m(1) agrees with
the Samaniego decomposition with the committed signature.

### Run on the cluster

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.precompile(); Pkg.test()'
    sbatch --array=1-12,14-18,20,22 templates/run_system.sbatch                       # n <= 27
    sbatch --array=13,19,21 --mem=24G --time=12:00:00 templates/run_system.sbatch     # n = 30
    grep -l "FAILED" logs/thresholds_*.out
    ls output/thresholds/row_*[0-9].csv | wc -l                                       # 682 when complete
    srun --partition=compute-slim --ntasks=1 --mem=1G --time=00:10:00 \
        julia --project=. scripts/summarize_thresholds.jl                             # output/thresholds_summary.csv

### Computational effort

`scripts/benchmark_sweep.jl` times, per system and on one core, the enumeration
needed for the signature of the complete system, the sweep (Algorithm alg), and the
sweep that also computes the MTTF of every state under the 31 models, and checks the
per-state MTTF against the classes (`output/benchmark/sweep_<name>.csv`):

    sbatch --array=1-12,14-18,20,22 templates/benchmark_sweep.sbatch
    sbatch --array=13,19,21 --mem=24G --time=12:00:00 templates/benchmark_sweep.sbatch

## Benchmark policies (no delay)

`src/policies.jl`, for every grid row and the same 45 cost pairs: the best
r-out-of-n policy (repair when r components have failed or the system failed),
the optimal count-based policy, and the optimal state-based policy π* of the
paper ("Optimal policy" in Section sec:policies: all state-based policies that
repair at failed states and continue at 1), all by Dinkelbach's method. The
state-based pass visits the 2^n states by levels, with level sums over sub-states
obtained from the sons (as in Lemma lemma:Wk), two levels in memory: about
4 (w+1) binomial(n, w) Float64 values (~0.5 GB for n = 23, ~10 GB for n = 27).

Output `output/policies/row_<id>.csv` (`# policies_version`), one line per cost pair:
`ltmc_best_gamma`, `best_gamma` (best MTTF-threshold policy), `best_r`, `ltmc_best_r`,
`ltmc_r1`, `ltmc_gamma1`, `ltmc_corrective`, `ltmc_count`, `count_policy` (`R` repair,
`.` continue, `-` unreachable, for f = 0..n-1), `count_threshold` (whether it is an
r-out-of-n policy), `ltmc_state`, `state_mixed_levels` (levels f where π* repairs at
some but not all working states), `state_repairs` (repairs/working states per
f = 1..n-1). State-based for n <= 23 in the main run (NaN otherwise); for
n = 24, 26, 27 (systems 20, 14, 17) in `output/policies_large/`; not for n = 30.

Every row checks state <= count <= best r-out-of-n and state <= best
MTTF-threshold, and that the r = n policy (computed with the count chain) equals
γ -> 0+ (computed from the classes), and r = 1 equals γ = 1 when γ = 1 continues
only at 1 (all systems except 6D, where 15 other states have m(x) = m(1) because
they differ from 1 only in irrelevant components).

Validation: `test/reference/policies_reference.py` follows the paper literally,
with explicit sums over all sub-states: policies evaluated with sys:m rep and
sys:c rep, π* by the algorithm of the paper and, for the bridge, by brute force
over all 2^15 state-based policies; count-based optimum by brute force over all
sets of counts; 72 cases (systems 1, 3, 9, 6D). `policies_port.py` is a Float64
port of the level-sum design (agrees to 1.4e-15). `Pkg.test()` compares the Julia
code with the reference (relative tolerance 1e-10) and checks that the state pass,
forced to follow π^γ, reproduces the threshold curves.

### Run on the cluster

    git pull
    module purge; module load softwares; module load Julia/1.10.4-linux-x86_64
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=4G --time=00:30:00 \
        julia --project=. -e 'using Pkg; Pkg.precompile(); Pkg.test()'
    # pilot: one row of system 22 (n = 23, the largest state-based case of the main run)
    srun --partition=compute-slim --ntasks=1 --cpus-per-task=1 --mem=4G --time=01:00:00 \
        julia --project=. scripts/run_policies.jl 1 1 --rows=652:652 --verbose
    sbatch templates/policies.sbatch                       # all 682 rows, 40 tasks
    grep -l "FAILED" logs/policies_*.out
    ls output/policies | wc -l                             # 682
    srun --partition=compute-slim --ntasks=1 --mem=1G --time=00:10:00 \
        julia --project=. scripts/summarize_policies.jl    # output/policies_summary.csv

The output of the pilot row is kept (the array job skips it).

## Layout

- `src/` reusable code (systems, signatures, failure models, per-state sweep, threshold policies, grid, benchmark policies)
- `scripts/` entry points
- `templates/` Slurm job templates
- `data/systems/` system definitions; `data/signatures/` computed signatures; `data/grid.csv` experiment grid (all committed)
- `test/reference/` independent reference values used by the tests
- `output/`, `logs/` results and job logs (git-ignored)
