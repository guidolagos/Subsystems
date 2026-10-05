# Per-state counting vectors: Lemma lemma:Wk and Algorithm alg of the paper.
#
# For a state x with w = |x| working components, W^x_k (k = 0..w) is the number of
# sub-states of x with k fewer working components in which the system works. The
# sweep visits the states by number of working components, w = 0, 1, ..., n
# (Algorithm alg runs over the number of failures k = n - w, from n down to 0), and
# obtains W^x from the vectors of the w sons of x (Lemma lemma:Wk):
#
#     W^x_0 = Φ(x),   W^x_k = (1/k) Σ_{x' ∈ Sons(x)} W^{x'}_{k-1},   k = 1..w,
#
# keeping in memory only the vectors of two consecutive levels. Φ is evaluated
# exactly once per state (2^n evaluations).
#
# States of a level are enumerated in increasing integer order, which is the
# colexicographic order of their sets of working components; the position (rank)
# of a son is computed from the bit positions of x, so no 2^n index arrays are
# needed. Level w is stored as a (w+1) x binomial(n, w) matrix (one column per
# state) of UInt16 or UInt32, depending on the largest possible entry binomial(w, w÷2).
#
# The working states are grouped into classes of equal counting vector (W^x has
# length w+1, so a class also fixes w). Everything the experiments need from a
# state (its MTTF under any LFMO model, and its contribution to the cost of the
# MTTF-threshold policies, src/thresholds.jl) depends on the state only through its
# class, so the sweep returns the classes with their multiplicities. A state x
# defines a minimal path set iff W^x_0 = 1 and W^x_1 = 0 (all sons failed).

const CLASSES_VERSION = 1
const MAX_SWEEP_N = 34          # entries binomial(w, w÷2) must fit in UInt32

"""
    ClassLevel

The classes of working states with `w` working components: column `c` of `W`
(a `(w+1) x nclasses` matrix) is the counting vector (W_0, ..., W_w) of class `c`,
and `mult[c]` is the number of states in it.
"""
struct ClassLevel
    W::Matrix{UInt32}
    mult::Vector{Int64}
end

"""
    SweepResult

Output of `sweep`: `levels[w+1]` are the classes of working states with `w`
working components; `working[w+1]` the number of such states; `mps` the number of
minimal path sets; `phi_evals` the number of evaluations of Φ; `seconds` the
duration of the sweep; `level_bytes` the peak memory of the two level matrices;
`mttf_sums[i]` the sum over all states of the MTTF computed per state with the
i-th `deltas` table (benchmark mode only, otherwise empty).
"""
struct SweepResult
    n::Int
    levels::Vector{ClassLevel}
    working::Vector{Int64}
    mps::Int64
    phi_evals::Int64
    seconds::Float64
    level_bytes::Int64
    mttf_sums::Vector{Float64}
end

nclasses(res::SweepResult) = sum(length(l.mult) for l in res.levels)

"binomial(p, j) for 0 <= p, j <= n, as B[p+1, j+1] (0 when j > p)."
binomial_table(n::Integer) = [binomial(p, j) for p in 0:n, j in 0:n]

level_eltype(w::Integer) = binomial(w, w ÷ 2) <= typemax(UInt16) ? UInt16 : UInt32

"Next integer with the same number of set bits (Gosper's hack)."
@inline function next_combination(x::UInt64)
    c = x & (~x + one(UInt64))
    r = x + c
    return (((r ⊻ x) >> 2) ÷ c) | r
end

"""
    sweep(cs::CompiledSystem; record = nothing, deltas = Matrix{Float64}[]) -> SweepResult

Algorithm alg: counting vectors of all 2^n states, grouped into classes.

- `record`: a `Dict{UInt64,Vector{Int}}` that receives the counting vector of every
  state (tests on small systems only).
- `deltas`: tables `D[w, k+1]` (as `ModelTables.delta`); for each one, the MTTF of
  every working state is computed on the fly as Σ_k W^x_k / binomial(w,k) D[w,k+1]
  and summed into `mttf_sums` (benchmark of the per-state computation).
"""
function sweep(cs::CompiledSystem; record::Union{Nothing,Dict{UInt64,Vector{Int}}} = nothing,
               deltas::Vector{Matrix{Float64}} = Matrix{Float64}[])
    n = cs.n
    n <= MAX_SWEEP_N || error("sweep: n = $(n) > $(MAX_SWEEP_N)")
    for D in deltas
        size(D, 1) >= n && size(D, 2) >= n || error("sweep: deltas table smaller than n x n")
    end
    B = binomial_table(n)
    t0 = time()
    works(cs, zero(UInt64)) && error("Φ(0,…,0) = 1: the system works with all components failed")
    record === nothing || (record[zero(UInt64)] = [0])
    prev = zeros(UInt16, 1, 1)                 # level 0: the failed state 0
    levels = ClassLevel[ClassLevel(zeros(UInt32, 1, 0), Int64[])]
    working = zeros(Int64, n + 1)
    mps = 0
    evals = 1
    sums = zeros(Float64, length(deltas))
    peak = sizeof(prev)
    for w in 1:n
        cur = zeros(level_eltype(w), w + 1, B[n + 1, w + 1])
        peak = max(peak, sizeof(prev) + sizeof(cur))
        lev, nw, nm = _sweep_level!(cur, prev, cs, w, B, record, deltas, sums)
        push!(levels, lev)
        working[w + 1] = nw
        mps += nm
        evals += size(cur, 2)
        prev = cur
    end
    seconds = time() - t0
    working[n + 1] == 1 || error("Φ(1,…,1) = 0: the system does not work with all components working")
    return SweepResult(n, levels, working, mps, evals, seconds, peak, sums)
end

function _sweep_level!(cur::Matrix{T}, prev::Matrix{S}, cs::CompiledSystem, w::Int,
                       B::Matrix{Int}, record, deltas::Vector{Matrix{Float64}},
                       sums::Vector{Float64}) where {T<:Unsigned,S<:Unsigned}
    nstates = size(cur, 2)
    pos = zeros(Int, w)              # bit positions (0-based) of the working components, increasing
    prefix = zeros(Int, w + 1)       # prefix[r] = Σ_{i<r} binomial(pos_i, i)
    suf = zeros(Int, w + 1)          # suf[r]    = Σ_{i>=r} binomial(pos_i, i-1)
    acc = zeros(Int64, w)
    index = Dict{Vector{T},Int}()
    reps = Vector{Vector{T}}()
    keybuf = zeros(T, w + 1)         # lookups use a Vector of the key type (same hash as the keys)
    mult = Int64[]
    nworking = 0
    nmps = 0
    x = (one(UInt64) << w) - one(UInt64)
    for idx in 1:nstates
        if works(cs, x)
            nworking += 1
            y = x
            i = 0
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
            fill!(acc, 0)
            for r in 1:w
                # colex rank (1-based column) of the son without the r-th working component
                c = prefix[r] + suf[r + 1] + 1
                @inbounds for j in 1:w
                    acc[j] += prev[j, c]
                end
            end
            cur[1, idx] = one(T)
            for j in 1:w
                q, rm = divrem(acc[j], j)
                rm == 0 || error("Lemma lemma:Wk: sum of the sons not divisible by $(j) (w = $(w))")
                @inbounds cur[j + 1, idx] = T(q)        # T(q) throws InexactError on overflow
            end
            cur[2, idx] == zero(T) && (nmps += 1)
            for j in 1:(w + 1)
                @inbounds keybuf[j] = cur[j, idx]
            end
            id = get(index, keybuf, 0)
            if id == 0
                key = copy(keybuf)
                push!(reps, key)
                push!(mult, 1)
                index[key] = length(reps)
            else
                mult[id] += 1
            end
            record === nothing || (record[x] = Int.(keybuf))
            for (mi, D) in enumerate(deltas)
                s = 0.0
                for k in 0:(w - 1)
                    @inbounds s += cur[k + 1, idx] / B[w + 1, k + 1] * D[w, k + 1]
                end
                sums[mi] += s
            end
        else
            record === nothing || (record[x] = zeros(Int, w + 1))
        end
        idx < nstates && (x = next_combination(x))
    end
    Wm = Matrix{UInt32}(undef, w + 1, length(reps))
    for c in eachindex(reps)
        Wm[:, c] .= reps[c]
    end
    return ClassLevel(Wm, mult), nworking, nmps
end

sweep(sys::TwoTerminalSystem; kwargs...) = sweep(compile_system(sys); kwargs...)

# ---------------------------------------------------------------------------
# Class files: output/classes/<name>.bin (binary, little-endian, git-ignored)
# ---------------------------------------------------------------------------

const CLASSES_MAGIC = 0x53554253595343c1    # arbitrary file signature

classes_path(sys::TwoTerminalSystem; dir::AbstractString = CLASSES_DIR) = joinpath(dir, "$(sys.name).bin")

"""
    write_classes(path, res::SweepResult)

Binary file: magic, version, n, mps, phi_evals, seconds, level_bytes, then for
w = 0..n: number of classes, number of working states, the UInt32 matrix W and the
Int64 multiplicities. Written atomically.
"""
function write_classes(path::AbstractString, res::SweepResult)
    io = IOBuffer()
    write(io, CLASSES_MAGIC, Int64(CLASSES_VERSION), Int64(res.n), Int64(res.mps),
          Int64(res.phi_evals), Float64(res.seconds), Int64(res.level_bytes))
    for w in 0:res.n
        lev = res.levels[w + 1]
        write(io, Int64(length(lev.mult)), Int64(res.working[w + 1]))
        write(io, lev.W)
        write(io, lev.mult)
    end
    return atomic_write(path, take!(io))
end

"""
    read_classes(path) -> SweepResult

Reads a file written by `write_classes`; refuses other versions. `mttf_sums` is empty.
"""
function read_classes(path::AbstractString)
    open(path, "r") do io
        read(io, UInt64) == CLASSES_MAGIC || error("$(path): not a class file")
        v = read(io, Int64)
        v == CLASSES_VERSION || error("$(path): class file version $(v), expected $(CLASSES_VERSION); recompute it")
        n = read(io, Int64); mps = read(io, Int64); evals = read(io, Int64)
        seconds = read(io, Float64); bytes = read(io, Int64)
        levels = ClassLevel[]
        working = zeros(Int64, n + 1)
        for w in 0:n
            nc = read(io, Int64)
            working[w + 1] = read(io, Int64)
            W = Matrix{UInt32}(undef, w + 1, nc)
            read!(io, W)
            mult = Vector{Int64}(undef, nc)
            read!(io, mult)
            sum(mult; init = 0) == working[w + 1] || error("$(path): multiplicities of level $(w) do not add up")
            push!(levels, ClassLevel(W, mult))
        end
        eof(io) || error("$(path): trailing data")
        return SweepResult(n, levels, working, mps, evals, seconds, bytes, Float64[])
    end
end

"""
    check_sweep(res, W) -> Bool

Consistency with the structural signature counts `W` of the complete system
(`W[k+1]` = working states with k failed components): the counting vector of the
state with all components working is `W`, and level w has `W[n-w+1]` working states.
"""
function check_sweep(res::SweepResult, W::AbstractVector{<:Integer})
    n = res.n
    length(W) == n + 1 || error("check_sweep: W has length $(length(W)), expected $(n + 1)")
    top = res.levels[n + 1]
    length(top.mult) == 1 && top.mult[1] == 1 || error("check_sweep: the top level must have one state")
    Int.(top.W[:, 1]) == collect(W) || error("check_sweep: counting vector of (1,…,1) differs from the signature counts")
    all(res.working[w + 1] == W[n - w + 1] for w in 0:n) || error("check_sweep: working states per level differ from W")
    return true
end
