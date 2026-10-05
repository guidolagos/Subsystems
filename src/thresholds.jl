# MTTF at every state and MTTF-threshold repair policies (Section sec:policies).
#
# Failure models are LFMO (src/models.jl). For a system with w components, let
# p^w_k be the probability that its survivor count visits w - k (k failed
# components); the count leaves that value at rate Ψ(w-k), so
#
#     δ^w_k := E T_{k+1:w} - E T_{k:w} = p^w_k / Ψ(w-k),   k = 0..w-1,
#
# and, by eq:E Samaniego sub, the MTTF at a state x with w working components is
#
#     m(x) = Σ_{k=0}^{w-1} W^x_k / binomial(w,k) δ^w_k.
#
# Both p and δ are computed with nonnegative terms only (no alternating sums).
#
# Policy π^γ repairs at the states with m(x) < γ m(1) and continues at the others.
# Under LFMO, m is monotone (x' <= x implies m(x') <= m(x): couple the two
# subsystems with the same subordinator and thresholds), so the continuation set
# C_γ = {x : m(x) >= γ m(1)} is an up-set. Every path from (1,…,1) to a state of
# C_γ stays in C_γ, so a state x ∈ C_γ is visited exactly as in the chain without
# repairs: its expected occupation time is τ_w = p^n_{n-w} / (binomial(n,w) Ψ(w)),
# independently of γ. Over one cycle (renewal at each repair):
#
#     E[cycle length]                 T = Σ_{x∈C_γ} τ_|x|
#     E[failed components at repair]  N = Σ_{x∈C_γ} τ_|x| |x| Ψ(1)
#     P(the repair is corrective)     P = Σ_{x∈C_γ} τ_|x| Σ_j q_|x|(j) (binomial(|x|,j) - W^x_j)
#
# where q_w(j) is the rate from a state with w working components to one given
# sub-state with j fewer (eq. def:lambda of the companion paper), so the last sum is
# the rate from x to failed states. With the repair cost
# c_set + (n - |y|) c_cmp + c_sys 1{Φ(y)=0} paid once per cycle,
#
#     LTMC^γ = (c_set + c_cmp N + c_sys P) / T.
#
# Each term depends on x only through its class (w, W^x) (src/sweep.jl), so a
# model is evaluated on the classes, sorted by m. The breakpoints are the distinct
# values m/m(1) > 0 (ties within RTOL_TIES m(1) merged); at breakpoint i the policy
# continues at the states with m(x) >= γ_i m(1), and LTMC^γ is constant for γ in
# (γ_{i+1}, γ_i]. The last breakpoint (γ -> 0+) is the corrective policy, for
# which T = m(1) and P = 1 (checked).
#
# Validation: test/reference/thresholds_reference.py solves the linear systems
# sys:m rep and sys:c rep of the paper directly, with m from the recursion eq:mttf
# (no signatures), and Pkg.test() compares (relative tolerance 1e-10).

const RTOL_TIES = 1e-12

"""
    ModelTables

Quantities of an LFMO model for systems of up to `n` components: `psi[w]` = Ψ(w);
`q[w, j]` the rate from a state with w working components to one given sub-state
with j fewer; `delta[w, k+1]` = δ^w_k = E T_{k+1:w} - E T_{k:w}; `tau[w]` the
expected time spent in one given state with w working components, starting from
all components working, without repairs.
"""
struct ModelTables
    n::Int
    psi::Vector{Float64}
    q::Matrix{Float64}
    delta::Matrix{Float64}
    tau::Vector{Float64}
end

function model_tables(m::LFMOModel, n::Integer)
    R, ex = transition_rates(m, n)          # indexed by failed components c = n - w
    psi = [ex[n - w + 1] for w in 1:n]
    q = zeros(n, n)
    P = zeros(n, n)                         # P[w, j]: j components fail at the next epoch
    for w in 1:n, j in 1:w
        r = R[n - w + 1, n - w + j + 1]
        P[w, j] = r / psi[w]
        q[w, j] = r / binomial(w, j)
    end
    delta = zeros(n, n)
    p = zeros(n)
    for w in 1:n
        fill!(p, 0.0)
        p[1] = 1.0                          # p[k+1] = P(the count of a w-system visits k failed)
        for k in 1:(w - 1)
            s = 0.0
            for i in 0:(k - 1)
                s += p[i + 1] * P[w - i, k - i]
            end
            p[k + 1] = s
        end
        for k in 0:(w - 1)
            delta[w, k + 1] = p[k + 1] / psi[w - k]
        end
    end
    tau = [delta[n, n - w + 1] / binomial(n, w) for w in 1:n]
    return ModelTables(n, psi, q, delta, tau)
end

"""
    expected_order_statistics(mt, w) -> Vector{Float64}

E T_{k:w}, k = 1..w, for a system of w components: partial sums of δ^w.
"""
expected_order_statistics(mt::ModelTables, w::Integer) = cumsum(mt.delta[w, 1:w])

"""
    class_mttf(mt, W, w) -> Float64

MTTF of a state with w working components and counting vector W (W[k+1] = W_k).
"""
function class_mttf(mt::ModelTables, W::AbstractVector, w::Integer)
    s = 0.0
    for k in 0:(w - 1)
        s += W[k + 1] / binomial(w, k) * mt.delta[w, k + 1]
    end
    return s
end

"""
    ThresholdCurve

Breakpoints of γ ↦ LTMC^γ for one system and model. `gamma` is decreasing,
`gamma[1] = 1`; `T`, `N`, `P` and `states` are cumulative over the classes with
m >= gamma[i] m(1): mean cycle length, mean number of failed components at the
repair, probability that the repair is corrective, and number of continuation states.
"""
struct ThresholdCurve
    n::Int
    m1::Float64
    gamma::Vector{Float64}
    T::Vector{Float64}
    N::Vector{Float64}
    P::Vector{Float64}
    states::Vector{Int64}
end

function threshold_curve(res::SweepResult, mt::ModelTables)
    n = res.n
    n <= mt.n || error("threshold_curve: model tables for n = $(mt.n) < $(n)")
    Bf = [Float64(binomial(w, j)) for w in 0:n, j in 0:n]
    total = sum(length(res.levels[w + 1].mult) for w in 1:n)
    m = Vector{Float64}(undef, total)
    t = similar(m); nf = similar(m); pf = similar(m)
    cnt = Vector{Int64}(undef, total)
    i = 0
    for w in 1:n
        lev = res.levels[w + 1]
        for c in eachindex(lev.mult)
            i += 1
            k = lev.mult[c]
            s = 0.0
            r = 0.0
            for kk in 0:(w - 1)
                s += lev.W[kk + 1, c] / Bf[w + 1, kk + 1] * mt.delta[w, kk + 1]
            end
            for j in 1:w
                r += mt.q[w, j] * (Bf[w + 1, j + 1] - lev.W[j + 1, c])
            end
            m[i] = s
            t[i] = k * mt.tau[w]
            nf[i] = k * mt.tau[w] * w * mt.psi[1]
            pf[i] = k * mt.tau[w] * r
            cnt[i] = k
        end
    end
    m1 = class_mttf(mt, view(res.levels[n + 1].W, :, 1), n)
    perm = sortperm(m; rev = true)
    m[perm[1]] <= m1 * (1 + RTOL_TIES) || error("threshold_curve: a state has larger MTTF than (1,…,1)")
    m[perm[end]] > 0 || error("threshold_curve: a working state has MTTF 0")
    gamma = Float64[]; Tc = Float64[]; Nc = Float64[]; Pc = Float64[]; Sc = Int64[]
    T = 0.0; N = 0.0; P = 0.0; S = 0
    rep = m[perm[1]]
    for i in perm
        if rep - m[i] > RTOL_TIES * m1          # new group: close the previous one
            push!(gamma, rep / m1); push!(Tc, T); push!(Nc, N); push!(Pc, P); push!(Sc, S)
            rep = m[i]
        end
        T += t[i]; N += nf[i]; P += pf[i]; S += cnt[i]
    end
    push!(gamma, rep / m1); push!(Tc, T); push!(Nc, N); push!(Pc, P); push!(Sc, S)
    abs(gamma[1] - 1) <= RTOL_TIES * 10 || error("threshold_curve: first breakpoint is not γ = 1")
    gamma[1] = 1.0
    # corrective policy (γ -> 0+): mean cycle = MTTF of the new system, repair always corrective
    isapprox(Tc[end], m1; rtol = 1e-9) || error("threshold_curve: T(0+) = $(Tc[end]) differs from m(1) = $(m1)")
    isapprox(Pc[end], 1.0; rtol = 1e-9) || error("threshold_curve: P(0+) = $(Pc[end]) differs from 1")
    Sc[end] == sum(res.working) || error("threshold_curve: the corrective policy must continue at every working state")
    return ThresholdCurve(n, m1, gamma, Tc, Nc, Pc, Sc)
end

"""
    curve_ltmc(cv; ratio, setup, c_sys = 1.0) -> Vector{Float64}

LTMC at each breakpoint, with c_cmp = ratio c_sys / n and c_set = setup c_sys.
"""
curve_ltmc(cv::ThresholdCurve; ratio::Real, setup::Real, c_sys::Real = 1.0) =
    (setup * c_sys .+ (ratio * c_sys / cv.n) .* cv.N .+ c_sys .* cv.P) ./ cv.T

"""
    local_minima(L; rtol = 1e-12) -> Int

Number of strict local minima of the sequence L after merging consecutive values
that are equal up to `rtol` (1 for a unimodal curve).
"""
function local_minima(L::AbstractVector{<:Real}; rtol::Real = 1e-12)
    v = Float64[]
    for x in L
        (isempty(v) || abs(x - v[end]) > rtol * max(abs(x), abs(v[end]))) && push!(v, x)
    end
    k = length(v)
    return count(i -> (i == 1 || v[i] < v[i - 1]) && (i == k || v[i] < v[i + 1]), 1:k)
end

"""
    curve_at(cv, γ) -> index

Index of the breakpoint that defines π^γ for an arbitrary γ in [0, 1]: the last
breakpoint with gamma >= γ (up to RTOL_TIES); γ = 0 means γ -> 0+ (corrective).
"""
function curve_at(cv::ThresholdCurve, γ::Real)
    γ <= 0 && return length(cv.gamma)
    i = findlast(g -> g >= γ - RTOL_TIES, cv.gamma)
    i === nothing && error("curve_at: γ = $(γ) > 1")
    return i
end
