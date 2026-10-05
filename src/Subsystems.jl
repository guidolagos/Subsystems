module Subsystems

# Paths of the data shipped with the repository, and of the (git-ignored) outputs.
const REPO_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(REPO_DIR, "data")
const SYSTEMS_DIR = joinpath(DATA_DIR, "systems")
const SIGNATURES_DIR = joinpath(DATA_DIR, "signatures")
const OUTPUT_DIR = joinpath(REPO_DIR, "output")
const CLASSES_DIR = joinpath(OUTPUT_DIR, "classes")
const SWEEP_DIR = joinpath(OUTPUT_DIR, "sweep")
const THRESH_DIR = joinpath(OUTPUT_DIR, "thresholds")
const BENCH_DIR = joinpath(OUTPUT_DIR, "benchmark")

include("systems.jl")
include("signatures.jl")
include("models.jl")
include("sweep.jl")
include("thresholds.jl")
include("grid.jl")

# systems
export TwoTerminalSystem, ncomponents, nnodes, load_systems, get_system, parse_yeh_raw
# structural signatures of the complete systems (exhaustive enumeration)
export CompiledSystem, compile_system, works, state_from_working,
       count_working_states, tail_signature, signature, check_counts,
       signature_path, write_signature_file, read_signature_file
# LFMO failure models calibrated by the pairwise beta factor
export MODELS, LFMOModel, lfmo_model, psi, transition_rates, mean_batch
# per-state counting vectors (Algorithm alg) and their classes
export ClassLevel, SweepResult, sweep, nclasses, classes_path, write_classes, read_classes,
       check_sweep, CLASSES_VERSION
# MTTF at all states and MTTF-threshold policies
export ModelTables, model_tables, class_mttf, expected_order_statistics, ThresholdCurve,
       threshold_curve, curve_ltmc, local_minima, curve_at
# experiment grid and runs
export GRID_PATH, GRID_BETAS, ROWS_PER_SYSTEM, COST_RATIOS, COST_SETUPS, CURVE_GAMMAS,
       THRESH_VERSION, THRESH_HEADER, make_grid, write_grid, read_grid, thresh_path, curve_path,
       breaks_path, read_threshold_output, write_threshold_outputs, sweep_system, run_system

end # module
