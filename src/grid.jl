# Experiment grid and per-system runs of the MTTF-threshold policies.
#
# data/grid.csv has one row per (system, failure model, beta): the 22 systems ×
# 31 models (Exp, and the six other LFMO models of src/models.jl with
# beta ∈ GRID_BETAS) = 682 rows; system k owns rows 31(k-1)+1 .. 31k.
#
# `run_system` sweeps the states of one system once (Algorithm alg, src/sweep.jl;
# the classes are cached in output/classes/<name>.bin) and then writes, for each of
# its rows, output/thresholds/:
#   row_<id>.csv         one line per cost pair (ratio, setup): best γ and LTMC;
#   row_<id>_curve.csv   T, N, P on the γ grid CURVE_GAMMAS (any cost pair can be
#                        evaluated from it);
#   row_<id>_breaks.csv  all breakpoints, when there are at most MAX_BREAKS_FILE.
# The main file is written last, atomically; a row whose main file exists is skipped.

const GRID_PATH = joinpath(DATA_DIR, "grid.csv")
const THRESH_VERSION = 1
const GRID_BETAS = ["0.01", "0.05", "0.1", "0.2", "0.5"]
const ROWS_PER_SYSTEM = 1 + (length(MODELS) - 1) * length(GRID_BETAS)
# cost pairs, as in the companion paper: ratio = n c_cmp / c_sys, setup = c_set / c_sys
const COST_RATIOS = [1 / 10, 1 / 5, 1 / 2, 2 / 3, 1.0, 4 / 3, 2.0, 4.0, 10.0]
const COST_SETUPS = [0.0, 0.05, 0.1, 0.2, 0.5]
const CURVE_GAMMAS = collect(0:200) ./ 200
const MAX_BREAKS_FILE = 10_000

const GRID_HEADER = ["row_id", "system", "model", "beta"]
const THRESH_HEADER = ["ratio", "setup", "n_breakpoints", "best_index", "best_gamma", "best_gamma_lo",
                       "ltmc_best", "ltmc_gamma1", "ltmc_corrective", "local_minima"]

"""
    make_grid(; family = "yeh2024") -> Vector{NamedTuple}

All (system, model, beta) rows, ordered by system, then model (`MODELS` order),
then beta.
"""
function make_grid(; family::AbstractString = "yeh2024")
    systems = load_systems(family)
    rows = NamedTuple{(:row_id, :system, :model, :beta),Tuple{Int,Int,String,String}}[]
    for sys in systems
        push!(rows, (row_id = length(rows) + 1, system = sys.id, model = "Exp", beta = "0"))
        for name in MODELS[2:end], b in GRID_BETAS
            push!(rows, (row_id = length(rows) + 1, system = sys.id, model = name, beta = b))
        end
    end
    return rows
end

function write_grid(path::AbstractString, rows)
    io = IOBuffer()
    println(io, "# Experiment grid of the Subsystems repository (scripts/make_grid.jl).")
    println(io, join(GRID_HEADER, ","))
    for r in rows
        println(io, r.row_id, ",", r.system, ",", r.model, ",", r.beta)
    end
    return atomic_write(path, take!(io))
end

function read_grid(path::AbstractString = GRID_PATH)
    header, rows = read_simple_csv(path)
    header == GRID_HEADER || error("$(path): unexpected header $(header)")
    out = [(row_id = parse(Int, r[1]), system = parse(Int, r[2]), model = r[3], beta = r[4]) for r in rows]
    [r.row_id for r in out] == collect(1:length(out)) || error("$(path): row ids must be 1..N")
    return out
end

thresh_path(id::Integer; dir::AbstractString = THRESH_DIR) = joinpath(dir, "row_$(id).csv")
curve_path(id::Integer; dir::AbstractString = THRESH_DIR) = joinpath(dir, "row_$(id)_curve.csv")
breaks_path(id::Integer; dir::AbstractString = THRESH_DIR) = joinpath(dir, "row_$(id)_breaks.csv")

_fmt(x::Real) = repr(Float64(x))

function _row_header(io, row, sys, cv, nclass, seconds)
    println(io, "# thresholds_version = ", THRESH_VERSION)
    println(io, "# row_id = ", row.row_id)
    println(io, "# system = ", sys.id)
    println(io, "# name = ", sys.name)
    println(io, "# model = ", row.model)
    println(io, "# beta = ", row.beta)
    println(io, "# n = ", cv.n)
    println(io, "# mttf_new = ", _fmt(cv.m1))
    println(io, "# n_breakpoints = ", length(cv.gamma))
    println(io, "# n_classes = ", nclass)
    println(io, "# seconds = ", round(seconds; digits = 4))
    println(io, "# julia_version = ", VERSION)
    println(io, "# host = ", _hostname())
end

"""
    write_threshold_outputs(row, sys, cv, nclass, seconds; dir = THRESH_DIR)

Writes the curve file, the breakpoints file (if small) and, last, the main file.
"""
function write_threshold_outputs(row, sys::TwoTerminalSystem, cv::ThresholdCurve, nclass::Integer,
                                 seconds::Real; dir::AbstractString = THRESH_DIR)
    io = IOBuffer()
    _row_header(io, row, sys, cv, nclass, seconds)
    println(io, "gamma,index,T,N,P,states")
    for g in CURVE_GAMMAS
        i = curve_at(cv, g)
        println(io, _fmt(g), ",", i, ",", _fmt(cv.T[i]), ",", _fmt(cv.N[i]), ",", _fmt(cv.P[i]), ",", cv.states[i])
    end
    atomic_write(curve_path(row.row_id; dir = dir), take!(io))
    if length(cv.gamma) <= MAX_BREAKS_FILE
        io = IOBuffer()
        _row_header(io, row, sys, cv, nclass, seconds)
        println(io, "index,gamma,T,N,P,states")
        for i in eachindex(cv.gamma)
            println(io, i, ",", _fmt(cv.gamma[i]), ",", _fmt(cv.T[i]), ",", _fmt(cv.N[i]), ",",
                    _fmt(cv.P[i]), ",", cv.states[i])
        end
        atomic_write(breaks_path(row.row_id; dir = dir), take!(io))
    end
    io = IOBuffer()
    _row_header(io, row, sys, cv, nclass, seconds)
    println(io, join(THRESH_HEADER, ","))
    for ratio in COST_RATIOS, setup in COST_SETUPS
        L = curve_ltmc(cv; ratio = ratio, setup = setup)
        ib = argmin(L)                       # first minimum: the largest γ among ties
        lo = ib < length(L) ? cv.gamma[ib + 1] : 0.0
        println(io, _fmt(ratio), ",", _fmt(setup), ",", length(L), ",", ib, ",", _fmt(cv.gamma[ib]), ",",
                _fmt(lo), ",", _fmt(L[ib]), ",", _fmt(L[1]), ",", _fmt(L[end]), ",", local_minima(L))
    end
    return atomic_write(thresh_path(row.row_id; dir = dir), take!(io))
end

"""
    read_threshold_output(path) -> (meta, header, rows)

Reads a main or curve file; refuses other versions.
"""
function read_threshold_output(path::AbstractString)
    meta = Dict{String,String}()
    for line in eachline(path)
        startswith(line, "#") || break
        m = match(r"^#\s*([A-Za-z_0-9]+)\s*=\s*(.*)$", line)
        m === nothing || (meta[m.captures[1]] = strip(m.captures[2]))
    end
    v = get(meta, "thresholds_version", "missing")
    v == string(THRESH_VERSION) || error("$(path): thresholds_version = $(v), expected $(THRESH_VERSION)")
    header, rows = read_simple_csv(path)
    return meta, header, rows
end

"""
    sweep_system(sys; classes_dir = CLASSES_DIR, sweep_dir = SWEEP_DIR, verbose = true) -> SweepResult

Reads the cached classes of `sys`, or runs the sweep, checks it against the
committed signature file, and caches it together with a one-line summary in
`<sweep_dir>/<name>.csv`.
"""
function sweep_system(sys::TwoTerminalSystem; classes_dir::AbstractString = CLASSES_DIR,
                      sweep_dir::AbstractString = SWEEP_DIR, verbose::Bool = true)
    path = classes_path(sys; dir = classes_dir)
    if isfile(path)
        res = read_classes(path)
        verbose && println("[classes] $(sys.name): read $(path) ($(nclasses(res)) classes)")
    else
        verbose && println("[sweep  ] $(sys.name): n = $(ncomponents(sys)), 2^n = $(2^ncomponents(sys)) states")
        flush(stdout)
        res = sweep(sys)
        write_classes(path, res)
        io = IOBuffer()
        println(io, "name,n,seconds,phi_evals,working_states,classes,mps,level_bytes,julia_version,host")
        println(io, sys.name, ",", res.n, ",", round(res.seconds; digits = 3), ",", res.phi_evals, ",",
                sum(res.working), ",", nclasses(res), ",", res.mps, ",", res.level_bytes, ",", VERSION, ",", _hostname())
        atomic_write(joinpath(sweep_dir, "$(sys.name).csv"), take!(io))
        verbose && println("[sweep  ] $(sys.name): $(round(res.seconds; digits = 2)) s, $(nclasses(res)) classes, ",
                           "$(res.mps) minimal path sets; wrote $(path)")
    end
    sig = signature_path(sys)
    if isfile(sig)
        check_sweep(res, read_signature_file(sig).W)
    else
        @warn "no signature file for $(sys.name); sweep not checked against it"
    end
    flush(stdout)
    return res
end

"""
    run_system(sysid; grid, systems, dir = THRESH_DIR, classes_dir, sweep_dir, rows = nothing) -> (done, skipped, failed)

All pending rows of one system (optionally only the row ids in `rows`).
"""
function run_system(sysid::Integer; grid = read_grid(), systems = load_systems(),
                    dir::AbstractString = THRESH_DIR, classes_dir::AbstractString = CLASSES_DIR,
                    sweep_dir::AbstractString = SWEEP_DIR, rows = nothing, verbose::Bool = true)
    sys = get_system(systems, sysid)
    mine = [r for r in grid if r.system == sysid && (rows === nothing || r.row_id in rows)]
    pending = [r for r in mine if !isfile(thresh_path(r.row_id; dir = dir))]
    skipped = length(mine) - length(pending)
    isempty(pending) && return (0, skipped, 0)
    res = sweep_system(sys; classes_dir = classes_dir, sweep_dir = sweep_dir, verbose = verbose)
    n = res.n
    s = Float64.(read_signature_file(signature_path(sys)).s)
    done = 0; failed = 0
    for r in pending
        try
            t0 = time()
            mt = model_tables(lfmo_model(r.model, r.beta), n)
            cv = threshold_curve(res, mt)
            # Samaniego: m(1) = Σ_k s_k E T_{k:n}
            et = expected_order_statistics(mt, n)
            isapprox(cv.m1, sum(s .* et); rtol = 1e-10) || error("m(1) differs from the Samaniego decomposition")
            seconds = time() - t0
            write_threshold_outputs(r, sys, cv, nclasses(res), seconds; dir = dir)
            done += 1
            verbose && println("[row $(r.row_id)] $(sys.name) $(r.model) β=$(r.beta): m(1) = ",
                               round(cv.m1; digits = 6), ", $(length(cv.gamma)) breakpoints, ",
                               round(seconds; digits = 3), " s")
        catch err
            failed += 1
            println("[row $(r.row_id)] FAILED: ", sprint(showerror, err))
        end
        flush(stdout)
    end
    return (done, skipped, failed)
end
