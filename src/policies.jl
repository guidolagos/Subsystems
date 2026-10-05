# Benchmarks for the MTTF-threshold policies (Section sec:policies, "Optimal
# policy"), with the repair starting as soon as it is decided (no delay):
#
#   * r-out-of-n policies: repair when r components have failed or the system
#     failed, r = 1..n;
#   * the optimal count-based policy: the action depends on the number f of
#     failed components (and the system working);
#   * the optimal state-based policy π*: over all state-based policies that repair
#     at every failed state and continue at 1.
#
# Every policy repairs at failed states. The cost of a repair with f failed
# components is c_set + f c_cmp, plus c_sys if the system failed (K(x) of the
# paper plus the setup c_set of the companion paper); repairing at 1 (f = 0) is not
# allowed. The long-term mean cost of a policy is E[cost of a cycle]/E[cycle
# length], and the optima are found by Dinkelbach's method (the algorithm of the
# paper): for a given g, a backward pass computes the policy minimizing
# E[C - g T] together with its E C and E T; then g := E C / E T, until
# E C - g E T >= -rtol max(1, |E C|).
#
# The structure is the one of src/optimal.jl of the companion repository
# Simple-delayed (no-delay case), re-written for the level sweep of src/sweep.jl:
#
#   count-based pass: backward recursion over f = n-1..0 with the survivor-count
#     chain and the tail signature Sbar_f (from the signature counts W);
#   state-based pass: the states are visited by number of working components
#     w = 0..n, in colex order within a level; for each state x the "level sums"
#     Σ_{y ⊆ x, |x|-|y| = j} C(y) (and the same for T) are obtained from those of
#     its sons divided by j (as in Lemma lemma:Wk, with C in place of Φ), keeping
#     two levels in memory; at failed states all sub-states are failed, so the
#     level sums are binomial(w,j) (c_set + (f+j) c_cmp + c_sys) and 0.
#
# Validation: test/reference/policies_reference.py (independent, follows the paper
# literally: explicit sub-state sums, policies evaluated by solving sys:m rep and
# sys:c rep, brute force over all state-based policies of the bridge and over all
# count-based policies); Pkg.test() compares. Each row also checks
# state <= count <= best r-out-of-n, state <= best MTTF-threshold policy, and
# that two pairs of policies computed by different code coincide: r = 1 with
# γ = 1 (when γ = 1 continues only at 1, i.e. no irrelevant components) and r = n
# with γ -> 0+ (corrective).

"""
    Costs(n; ratio, setup, c_sys = 1.0)

c_cmp = ratio c_sys / n and c_set = setup c_sys.
"""
struct Costs
    c_sys::Float64
    c_set::Float64
    c_cmp::Float64
end

Costs(n::Integer; ratio::Real, setup::Real, c_sys::Real = 1.0) =
    Costs(Float64(c_sys), Float64(setup * c_sys), Float64(ratio * c_sys / n))

repair_cost(c::Costs, f::Integer) = c.c_set + c.c_cmp * f
failure_cost(c::Costs, f::Integer) = c.c_set + c.c_cmp * f + c.c_sys

"""
    dinkelbach(pass, g0; maxit = 100, rtol = 1e-12) -> (g, policy, iterations)

`pass(g)` returns (policy, E C, E T) of a policy minimizing E[C - g T].
"""
function dinkelbach(pass, g0::Real; maxit::Integer = 100, rtol::Real = 1e-12)
    g = Float64(g0)
    for it in 1:maxit
        pol, c, l = pass(g)
        l > 0 || error("Dinkelbach: expected cycle length $(l) <= 0")
        c - g * l >= -rtol * max(1.0, abs(c)) && return (c / l, pol, it)
        g = c / l
    end
    error("Dinkelbach did not converge in $(maxit) iterations")
end

# ---------------------------------------------------------------------------
# Count-based policies
# ---------------------------------------------------------------------------

"""
    CountModel(W, mt)

Tail signature `Sbar[f+1] = W_f / binomial(n, f)` of the system and the model
tables `mt` (src/thresholds.jl) for n components.
"""
struct CountModel
    n::Int
    Sbar::Vector{Float64}
    P::Matrix{Float64}          # P[w, j]: probability that j of w working components fail at the next epoch
    psi::Vector{Float64}
end

function CountModel(W::AbstractVector{<:Integer}, mt::ModelTables)
    n = length(W) - 1
    mt.n >= n || error("CountModel: model tables for n = $(mt.n) < $(n)")
    Sbar = [Float64(big(W[f + 1]) // binomial(big(n), f)) for f in 0:n]
    P = [j <= w ? binomial(w, j) * mt.q[w, j] / mt.psi[w] : 0.0 for w in 1:n, j in 1:n]
    return CountModel(n, Sbar, P, mt.psi[1:n])
end

"""
    count_pass(cm, costs, g; forced = nothing) -> (actions, E C, E T)

Backward recursion over f = n-1..0. `actions[f+1]` is 1 (repair), 0 (continue)
or -1 (Sbar_f = 0, unreachable). With `forced = f -> Bool` the actions at f >= 1
are fixed, which evaluates that policy (e.g. `f -> f >= r`). Ties: repair.
"""
function count_pass(cm::CountModel, costs::Costs, g::Real; forced = nothing)
    n = cm.n; Sb = cm.Sbar
    c = zeros(n + 1); l = zeros(n + 1)
    act = fill(Int8(-1), n)
    for f in (n - 1):-1:0
        Sb[f + 1] > 0 || continue
        w = n - f
        c_cont = 0.0; l_cont = 1 / cm.psi[w]
        for j in 1:w
            Pj = cm.P[w, j]
            Pj == 0 && continue
            ps = Sb[f + j + 1] / Sb[f + 1]           # the system still works after j more failures
            c_cont += Pj * (ps * c[f + j + 1] + (1 - ps) * failure_cost(costs, f + j))
            l_cont += Pj * ps * l[f + j + 1]
        end
        call = if f == 0
            false
        elseif forced !== nothing
            forced(f)::Bool
        else
            repair_cost(costs, f) <= c_cont - g * l_cont
        end
        act[f + 1] = call ? Int8(1) : Int8(0)
        c[f + 1] = call ? repair_cost(costs, f) : c_cont
        l[f + 1] = call ? 0.0 : l_cont
    end
    return act, c[1], l[1]
end

"""
    rout_ltmc(cm, costs) -> Vector{Float64}

Long-term mean cost of the r-out-of-n policy, r = 1..n (repair when f >= r or the
system failed).
"""
function rout_ltmc(cm::CountModel, costs::Costs)
    return [begin
                _, c, l = count_pass(cm, costs, 0.0; forced = f -> f >= r)
                c / l
            end for r in 1:cm.n]
end

"Whether a count-based policy is an r-out-of-n policy on the reachable counts f >= 1."
function is_threshold(act::AbstractVector{Int8})
    seen = false
    for f in 1:(length(act) - 1)
        a = act[f + 1]
        a < 0 && continue
        if a == 1
            seen = true
        elseif seen
            return false
        end
    end
    return true
end

"String of a count policy: 'R' repair, '.' continue, '-' unreachable, for f = 0..n-1."
policy_string(act::AbstractVector{Int8}) = String([a < 0 ? '-' : (a == 1 ? 'R' : '.') for a in act])

# ---------------------------------------------------------------------------
# State-based policies
# ---------------------------------------------------------------------------

"""
    StatePhi(cs)

Structure function at all 2^n states (`phi[x+1]`), evaluated once per system.
"""
struct StatePhi
    n::Int
    phi::BitVector
end

function StatePhi(cs::CompiledSystem)
    n = cs.n
    n <= 30 || error("StatePhi: n = $(n) > 30")
    phi = falses(1 << n)
    for x in 0:((1 << n) - 1)
        phi[x + 1] = works(cs, UInt64(x))
    end
    return StatePhi(n, phi)
end

"""
    state_pass(sp, mt, costs, g; forced = nothing) -> ((repairs, working), E C, E T)

One pass of the algorithm of the paper ("Optimal policy"): for every state, in
increasing order of the number of working components, the action minimizing
E[C - g T] until the next repair. `repairs[f+1]` and `working[f+1]` count, among
the working states with f failed components, those where the policy repairs.
`forced = (x, f) -> Bool` fixes the actions at working states x ≠ 1 (evaluation of
a given policy). Ties: repair.
"""
function state_pass(sp::StatePhi, mt::ModelTables, costs::Costs, g::Real; forced = nothing)
    n = sp.n
    mt.n >= n || error("state_pass: model tables for n = $(mt.n) < $(n)")
    B = binomial_table(n)
    Bf = Float64.(B)
    repairs = zeros(Int, n + 1); working = zeros(Int, n + 1)
    pj = zeros(n)
    pos = zeros(Int, n); prefix = zeros(Int, n + 1); suf = zeros(Int, n + 1)
    # level 0: the state 0 is failed (Φ(0) = 0 for a semi-coherent system)
    sp.phi[1] && error("state_pass: Φ(0) = 1")
    prevC = fill(failure_cost(costs, n), 1, 1)
    prevL = zeros(1, 1)
    for w in 1:n
        f = n - w
        m = B[n + 1, w + 1]
        curC = Matrix{Float64}(undef, w + 1, m)
        curL = Matrix{Float64}(undef, w + 1, m)
        for j in 1:w
            pj[j] = mt.q[w, j] / mt.psi[w]           # to one given sub-state with j fewer working
        end
        x = (one(UInt64) << w) - one(UInt64)
        for idx in 1:m
            if !sp.phi[x + 1]
                for j in 0:w
                    @inbounds curC[j + 1, idx] = Bf[w + 1, j + 1] * failure_cost(costs, f + j)
                    @inbounds curL[j + 1, idx] = 0.0
                end
            else
                y = x; i = 0
                while y != zero(UInt64)
                    i += 1
                    @inbounds pos[i] = trailing_zeros(y)
                    y &= y - one(UInt64)
                end
                prefix[1] = 0
                for i in 1:w
                    @inbounds prefix[i + 1] = prefix[i] + B[pos[i] + 1, i + 1]
                end
                suf[w + 1] = 0
                for i in w:-1:1
                    @inbounds suf[i] = suf[i + 1] + B[pos[i] + 1, i]
                end
                for j in 1:w
                    @inbounds curC[j + 1, idx] = 0.0
                    @inbounds curL[j + 1, idx] = 0.0
                end
                for r in 1:w
                    col = prefix[r] + suf[r + 1] + 1      # colex rank of the son without the r-th working component
                    @inbounds for j in 1:w
                        curC[j + 1, idx] += prevC[j, col]
                        curL[j + 1, idx] += prevL[j, col]
                    end
                end
                sc = 0.0; sl = 0.0
                @inbounds for j in 1:w
                    curC[j + 1, idx] /= j
                    curL[j + 1, idx] /= j
                    sc += pj[j] * curC[j + 1, idx]
                    sl += pj[j] * curL[j + 1, idx]
                end
                c_cont = sc
                l_cont = 1 / mt.psi[w] + sl
                call = if f == 0
                    false
                elseif forced !== nothing
                    forced(x, f)::Bool
                else
                    repair_cost(costs, f) <= c_cont - g * l_cont
                end
                working[f + 1] += 1
                repairs[f + 1] += call
                curC[1, idx] = call ? repair_cost(costs, f) : c_cont
                curL[1, idx] = call ? 0.0 : l_cont
            end
            idx < m && (x = next_combination(x))
        end
        prevC = curC; prevL = curL
    end
    return (repairs, working), prevC[1, 1], prevL[1, 1]
end

"String of a state policy: repairs/working states for f = 1..n-1, separated by ';'."
state_calls_string(repairs, working) = join(("$(repairs[f + 1])/$(working[f + 1])" for f in 1:(length(repairs) - 2)), ";")

"Number of levels f at which the state-based policy repairs at some, but not all, working states."
mixed_levels(repairs, working) = count(f -> 0 < repairs[f] < working[f], eachindex(repairs))

# ---------------------------------------------------------------------------
# One grid row
# ---------------------------------------------------------------------------

const POLICIES_VERSION = 1
const POLICIES_HEADER = ["ratio", "setup", "ltmc_best_gamma", "best_gamma", "best_r", "ltmc_best_r",
                         "ltmc_r1", "ltmc_gamma1", "ltmc_corrective", "ltmc_count", "iters_count",
                         "count_policy", "count_threshold", "ltmc_state", "iters_state",
                         "state_mixed_levels", "state_repairs", "seconds_state"]

policies_path(id::Integer; dir::AbstractString = POLICIES_DIR) = joinpath(dir, "row_$(id).csv")

"""
    run_policies_row(row; systems, dir = POLICIES_DIR, state_max_n = 23, cache = Dict(), verbose = false) -> Symbol

For one grid row and the 45 cost pairs: best MTTF-threshold policy (from the
classes, src/thresholds.jl), r-out-of-n policies, optimal count-based policy and,
if n <= `state_max_n`, optimal state-based policy (Dinkelbach started at the best
MTTF-threshold policy). Writes `dir/row_<id>.csv` atomically; `:skipped` if it exists.
"""
function run_policies_row(row; systems::Vector{TwoTerminalSystem}, dir::AbstractString = POLICIES_DIR,
                          classes_dir::AbstractString = CLASSES_DIR, sweep_dir::AbstractString = SWEEP_DIR,
                          state_max_n::Integer = 23, cache::Dict = Dict(), verbose::Bool = false)
    path = policies_path(row.row_id; dir = dir)
    isfile(path) && return :skipped
    t0 = time()
    sys = get_system(systems, row.system)
    n = ncomponents(sys)
    sigfile = get!(() -> read_signature_file(signature_path(sys)), cache, (:sig, sys.id))
    res = get!(() -> sweep_system(sys; classes_dir = classes_dir, sweep_dir = sweep_dir, verbose = verbose),
               cache, (:classes, sys.id))
    mt = model_tables(lfmo_model(row.model, row.beta), n)
    cv = threshold_curve(res, mt)
    cm = CountModel(sigfile.W, mt)
    # γ = 1 continues only at 1 (so it is the r = 1 policy) unless some component is irrelevant
    gamma1_is_r1 = cv.states[1] == 1
    do_state = n <= state_max_n
    sp = do_state ? get!(() -> StatePhi(compile_system(sys)), cache, (:phi, sys.id)) : nothing
    if do_state
        # keep only the structure function of the current system in memory
        for k in collect(keys(cache))
            k isa Tuple && k[1] == :phi && k[2] != sys.id && delete!(cache, k)
        end
    end
    tol = 1e-9
    lines = String[]
    for ρ in COST_RATIOS, σ in COST_SETUPS
        costs = Costs(n; ratio = ρ, setup = σ)
        Lγ = curve_ltmc(cv; ratio = ρ, setup = σ)
        iγ = argmin(Lγ)
        Lr = rout_ltmc(cm, costs)
        ir = argmin(Lr)
        # two pairs of identical policies computed by different code
        isapprox(Lr[n], Lγ[end]; rtol = 1e-9) || error("row $(row.row_id): r = n ($(Lr[n])) differs from γ -> 0+ ($(Lγ[end]))")
        !gamma1_is_r1 || isapprox(Lr[1], Lγ[1]; rtol = 1e-9) ||
            error("row $(row.row_id): r = 1 ($(Lr[1])) differs from γ = 1 ($(Lγ[1]))")
        gc, actc, itc = dinkelbach(g -> count_pass(cm, costs, g), Lr[ir])
        gc <= Lr[ir] * (1 + tol) || error("row $(row.row_id): count-based optimum $(gc) > best r-out-of-n $(Lr[ir])")
        gs = NaN; its = 0; mixed = ""; reps = ""; secs = 0.0
        if do_state
            ts = time()
            gs, (rep, wk), its = dinkelbach(g -> state_pass(sp, mt, costs, g), min(Lγ[iγ], gc))
            secs = time() - ts
            gs <= gc * (1 + tol) || error("row $(row.row_id): state-based optimum $(gs) > count-based $(gc)")
            gs <= Lγ[iγ] * (1 + tol) || error("row $(row.row_id): state-based optimum $(gs) > best MTTF-threshold $(Lγ[iγ])")
            mixed = string(mixed_levels(rep, wk))
            reps = state_calls_string(rep, wk)
            verbose && println("  ratio $(ρ) setup $(σ): state-based $(its) passes, $(round(secs; digits = 2)) s")
        end
        push!(lines, join((repr(ρ), repr(σ), repr(Lγ[iγ]), repr(cv.gamma[iγ]), ir, repr(Lr[ir]), repr(Lr[1]),
                           repr(Lγ[1]), repr(Lγ[end]), repr(gc), itc, policy_string(actc), is_threshold(actc),
                           repr(gs), its, mixed, reps, round(secs; digits = 3)), ","))
    end
    seconds = time() - t0
    io = IOBuffer()
    println(io, "# policies_version = ", POLICIES_VERSION)
    println(io, "# row_id = ", row.row_id)
    println(io, "# system = ", sys.id)
    println(io, "# name = ", sys.name)
    println(io, "# model = ", row.model)
    println(io, "# beta = ", row.beta)
    println(io, "# n = ", n)
    println(io, "# states_gamma1 = ", cv.states[1])
    println(io, "# state_based = ", do_state)
    println(io, "# seconds = ", round(seconds; digits = 3))
    println(io, "# julia_version = ", VERSION)
    println(io, "# host = ", _hostname())
    println(io, join(POLICIES_HEADER, ","))
    for line in lines
        println(io, line)
    end
    atomic_write(path, take!(io))
    return :done
end

"""
    read_policies_output(path) -> (meta, header, rows)
"""
function read_policies_output(path::AbstractString)
    meta = Dict{String,String}()
    for line in eachline(path)
        startswith(line, "#") || break
        m = match(r"^#\s*([A-Za-z_]+)\s*=\s*(.*)$", line)
        m === nothing || (meta[m.captures[1]] = strip(m.captures[2]))
    end
    get(meta, "policies_version", "") == string(POLICIES_VERSION) ||
        error("$(path): policies_version $(get(meta, "policies_version", "?")), expected $(POLICIES_VERSION)")
    header, rows = read_simple_csv(path)
    header == POLICIES_HEADER || error("$(path): unexpected header")
    return meta, header, rows
end
