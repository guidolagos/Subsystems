module Subsystems

# Paths of the data shipped with the repository.
const REPO_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(REPO_DIR, "data")
const SYSTEMS_DIR = joinpath(DATA_DIR, "systems")
const SIGNATURES_DIR = joinpath(DATA_DIR, "signatures")

include("systems.jl")
include("signatures.jl")
include("models.jl")

# systems
export TwoTerminalSystem, ncomponents, nnodes, load_systems, get_system, parse_yeh_raw
# structural signatures of the complete systems (exhaustive enumeration)
export CompiledSystem, compile_system, works, state_from_working,
       count_working_states, tail_signature, signature, check_counts,
       signature_path, write_signature_file, read_signature_file
# LFMO failure models calibrated by the pairwise beta factor
export MODELS, LFMOModel, lfmo_model, psi, transition_rates, mean_batch

end # module
