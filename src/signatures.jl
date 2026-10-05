# Structural signatures by exhaustive enumeration of the 2^n states.
#
# State convention: a state is an unsigned integer x with bit (e-1) equal to 1
# iff component e WORKS. For each state we decide whether the source can reach
# the target through working components, and we tally
#
#     W[k+1] = #{ x : system works in x, and x has exactly k failed components },
#
# for k = 0..n. Then (see e.g. Marichal et al. 2011, and eq. (def:s) of the paper)
#
#     tail signature  Sbar_k = W_k / binomial(n, k),   k = 0..n,
#     signature       s_k    = Sbar_{k-1} - Sbar_k,     k = 1..n.
#
# All counts are exact integers and the signature is returned as exact
# rationals (Rational{BigInt}).

const MAX_COMPONENTS = 40   # 2^40 states would take days; guard against typos
const MAX_NODES = 64        # node sets are stored as bits of a UInt64

"""
    CompiledSystem

Bitmask representation of a `TwoTerminalSystem` used in the inner loop:
`from_mask[e]` and `to_mask[e]` are the one-bit node masks of the endpoints
of component `e`.
"""
struct CompiledSystem
    n::Int
    directed::Bool
    from_mask::Vector{UInt64}
    to_mask::Vector{UInt64}
    source_mask::UInt64
    target_mask::UInt64
end

function compile_system(sys::TwoTerminalSystem)
    n = ncomponents(sys)
    1 <= n <= MAX_COMPONENTS || error("$(sys.name): n = $(n) components is outside 1..$(MAX_COMPONENTS)")
    index = Dict{String,Int}()          # node name => 0-based bit position
    for (u, v) in sys.edges
        haskey(index, u) || (index[u] = length(index))
        haskey(index, v) || (index[v] = length(index))
    end
    length(index) <= MAX_NODES || error("$(sys.name): more than $(MAX_NODES) nodes")
    haskey(index, sys.source) || error("$(sys.name): source $(sys.source) is not an endpoint of any component")
    haskey(index, sys.target) || error("$(sys.name): target $(sys.target) is not an endpoint of any component")
    sys.source != sys.target || error("$(sys.name): source and target coincide")
    bit(name) = UInt64(1) << index[name]
    return CompiledSystem(
        n,
        sys.directed,
        [bit(u) for (u, _) in sys.edges],
        [bit(v) for (_, v) in sys.edges],
        bit(sys.source),
        bit(sys.target),
    )
end

"""
    works(cs::CompiledSystem, x::UInt64) -> Bool

Structure function: true iff the target is reachable from the source using
only the components whose bit is set in `x`. Computes the reachable node set
by repeated sweeps over the working components until it stops growing.
"""
@inline function works(cs::CompiledSystem, x::UInt64)
    reach = cs.source_mask
    while true
        previous = reach
        y = x
        while y != zero(UInt64)
            e = trailing_zeros(y) + 1
            y &= y - one(UInt64)
            @inbounds u = cs.from_mask[e]
            @inbounds v = cs.to_mask[e]
            if reach & u != zero(UInt64)
                reach |= v
            end
            if !cs.directed && reach & v != zero(UInt64)
                reach |= u
            end
        end
        reach & cs.target_mask != zero(UInt64) && return true
        reach == previous && return false
    end
end

works(sys::TwoTerminalSystem, x::Integer) = works(compile_system(sys), UInt64(x))

"""
    state_from_working(n, working_components) -> UInt64

State in which exactly the listed components (1-based) work.
"""
function state_from_working(n::Integer, working)
    x = zero(UInt64)
    for e in working
        1 <= e <= n || error("component $(e) outside 1..$(n)")
        x |= UInt64(1) << (e - 1)
    end
    return x
end

function _count_range(cs::CompiledSystem, lo::UInt64, hi::UInt64)
    W = zeros(Int, cs.n + 1)
    x = lo
    while x < hi
        if works(cs, x)
            @inbounds W[cs.n - count_ones(x) + 1] += 1
        end
        x += one(UInt64)
    end
    return W
end

"""
    count_working_states(sys; nchunks = 8 * Threads.nthreads()) -> Vector{Int}

Exhaustive enumeration of the 2^n states. Returns `W` with `W[k+1]` the
number of working states with exactly `k` failed components, `k = 0..n`.
The state space is split into `nchunks` contiguous ranges processed as
parallel tasks; the result does not depend on the number of threads.
"""
function count_working_states(cs::CompiledSystem; nchunks::Integer = 8 * Threads.nthreads())
    total = UInt64(1) << cs.n
    m = Int(min(UInt64(max(nchunks, 1)), total))
    bounds = [div(total * UInt64(i), UInt64(m)) for i in 0:m]
    tasks = [Threads.@spawn(_count_range(cs, bounds[i], bounds[i + 1])) for i in 1:m]
    W = zeros(Int, cs.n + 1)
    for t in tasks
        W .+= fetch(t)::Vector{Int}
    end
    sum(W) <= total || error("internal error: more working states than states")
    return W
end

count_working_states(sys::TwoTerminalSystem; kwargs...) =
    count_working_states(compile_system(sys); kwargs...)

"""
    tail_signature(W) -> Vector{Rational{BigInt}}

`Sbar[k+1] = W[k+1] / binomial(n, k)`, for k = 0..n.
"""
function tail_signature(W::AbstractVector{<:Integer})
    n = length(W) - 1
    return [big(W[k + 1]) // binomial(big(n), k) for k in 0:n]
end

"""
    signature(W) -> Vector{Rational{BigInt}}

Structural signature `s[k] = Sbar_{k-1} - Sbar_k`, for k = 1..n.
Throws an error unless the system is semi-coherent (Φ(1)=1, Φ(0)=0 and
monotone, which makes every s[k] nonnegative).
"""
function signature(W::AbstractVector{<:Integer})
    check_counts(W)
    Sbar = tail_signature(W)
    n = length(W) - 1
    s = [Sbar[k] - Sbar[k + 1] for k in 1:n]
    all(>=(0), s) || error("negative signature entry: the structure function is not monotone")
    sum(s) == 1 || error("signature does not sum to one")
    return s
end

"""
    check_counts(W)

Sanity checks on the counts: `W[1] == 1` (system works when all components
work), `W[end] == 0` (system fails when all fail) and
`0 <= W[k+1] <= binomial(n, k)`.
"""
function check_counts(W::AbstractVector{<:Integer})
    n = length(W) - 1
    n >= 1 || error("empty count vector")
    W[1] == 1 || error("Φ(1,…,1) = 0: the system does not work with all components working")
    W[end] == 0 || error("Φ(0,…,0) = 1: the system works with all components failed")
    for k in 0:n
        0 <= W[k + 1] <= binomial(big(n), k) || error("W[$(k)] = $(W[k + 1]) is out of range")
    end
    return true
end

# ---------------------------------------------------------------------------
# Signature files
# ---------------------------------------------------------------------------

signature_path(sys::TwoTerminalSystem; dir::AbstractString = SIGNATURES_DIR) =
    joinpath(dir, "$(sys.name).csv")

ratstring(r::Rational) = string(numerator(r), "/", denominator(r))

function parse_rational(s::AbstractString)
    parts = split(s, '/')
    length(parts) == 2 || error("not a rational: $(s)")
    return parse(BigInt, parts[1]) // parse(BigInt, parts[2])
end

"""
    atomic_write(path, content)

Write `content` to a temporary file in the same directory and rename it to
`path`, so that `path` is either absent or complete.
"""
function atomic_write(path::AbstractString, content)
    mkpath(dirname(abspath(path)))
    tmp = string(path, ".tmp.", getpid())
    try
        write(tmp, content)
        mv(tmp, path; force = true)
    finally
        isfile(tmp) && rm(tmp; force = true)
    end
    return path
end

function _hostname()
    try
        return Libc.gethostname()
    catch
        return "unknown"
    end
end

"""
    write_signature_file(path, sys, W; seconds, threads)

Write the counts, tail signature and signature of `sys` to `path` (atomic).
Columns: `k, W, binomial, tail_signature, signature, signature_float`, one row
per k = 0..n; exact values are written as `numerator/denominator`, and
`signature` at k = 0 is 0 by convention.
"""
function write_signature_file(path::AbstractString, sys::TwoTerminalSystem, W::AbstractVector{<:Integer};
                              seconds::Real = NaN, threads::Integer = Threads.nthreads())
    n = ncomponents(sys)
    length(W) == n + 1 || error("count vector has wrong length")
    s = signature(W)                       # also runs check_counts
    Sbar = tail_signature(W)
    io = IOBuffer()
    println(io, "# Subsystems structural signature")
    println(io, "# system = ", sys.id)
    println(io, "# name = ", sys.name)
    println(io, "# label = ", sys.label)
    println(io, "# directed = ", sys.directed)
    println(io, "# source = ", sys.source)
    println(io, "# target = ", sys.target)
    println(io, "# n_components = ", n)
    println(io, "# n_nodes = ", nnodes(sys))
    println(io, "# method = exhaustive enumeration of all 2^n states (source-target connectivity)")
    println(io, "# seconds = ", round(seconds; digits = 3))
    println(io, "# threads = ", threads)
    println(io, "# julia_version = ", VERSION)
    println(io, "# host = ", _hostname())
    println(io, "k,W,binomial,tail_signature,signature,signature_float")
    for k in 0:n
        sk = k == 0 ? big(0) // big(1) : s[k]
        println(io, k, ",", W[k + 1], ",", binomial(big(n), k), ",",
                ratstring(Sbar[k + 1]), ",", ratstring(sk), ",", Float64(sk))
    end
    return atomic_write(path, take!(io))
end

"""
    read_signature_file(path) -> NamedTuple

Returns `(meta = Dict{String,String}, W = Vector{Int}, s = Vector{Rational{BigInt}})`,
with `s[k]` for k = 1..n. Checks internal consistency of the file.
"""
function read_signature_file(path::AbstractString)
    meta = Dict{String,String}()
    for line in eachline(path)
        startswith(line, "#") || continue
        m = match(r"^#\s*([A-Za-z_0-9]+)\s*=\s*(.*)$", line)
        m === nothing || (meta[m.captures[1]] = strip(m.captures[2]))
    end
    header, rows = read_simple_csv(path)
    ck, cW, cs = column(header, "k"), column(header, "W"), column(header, "signature")
    [parse(Int, r[ck]) for r in rows] == collect(0:length(rows) - 1) || error("$(path): bad k column")
    W = [parse(Int, r[cW]) for r in rows]
    s = [parse_rational(r[cs]) for r in rows[2:end]]
    s == signature(W) || error("$(path): signature column inconsistent with W column")
    return (meta = meta, W = W, s = s)
end
