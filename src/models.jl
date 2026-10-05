# LFMO failure models calibrated by the pairwise beta factor.
#
# Every model is an LFMO distribution with Laplace exponent Ψ normalized to
# Ψ(1) = 1 (unit mean component lifetime) and a single free parameter, chosen
# so that the pairwise beta factor
#
#     β := (2Ψ(1) - Ψ(2)) / Ψ(1) = 2 - Ψ(2)
#
# equals a target value. β is the fraction of a component's failure rate due to
# events that also kill one given other component (see the paper, Section
# "Failure models"). All Ψ values and transition rates are computed in BigFloat,
# because the rates are alternating binomial sums.
#
#   model      parameter  Ψ(x), x > 0                          β ↔ parameter
#   Exp        —          x                                    β = 0
#   LFMO-C     p          x / (1 + p(x-1))                     p = β/(2-β)
#   LFMO-G     η          log(1+x/η) / log(1+1/η)              bisection (β decreasing in η)
#   LFMO-IG    η          (√(2x+η²)-η) / (√(2+η²)-η)           bisection; needs β < 2-√2
#   LFMO-S     α          x^α                                  α = log2(2-β)
#   BFR        p          (1-(1-p)^x) / p                      p = β
#   LSh        θ          (1-θ)x + θ                           θ = β
#
# Ψ(0) := 0 for every model ([LSh] is a killed subordinator).

const MODELS = ("Exp", "LFMO-C", "LFMO-G", "LFMO-IG", "LFMO-S", "BFR", "LSh")
const BIG_PRECISION = 512   # bits

"""
    LFMOModel

A calibrated failure model: `name` (one of `MODELS`), target pairwise beta
factor `beta` (given as a decimal string, e.g. "0.05", and parsed exactly) and
the calibrated parameter `param` (BigFloat; 0 for Exp).
"""
struct LFMOModel
    name::String
    beta::String
    param::BigFloat
end

function _psi_big(name::AbstractString, a::BigFloat, x::Integer)
    x == 0 && return BigFloat(0)
    X = BigFloat(x)
    name == "Exp" && return X
    name == "LFMO-C" && return X / (1 + a * (X - 1))
    name == "LFMO-G" && return log1p(X / a) / log1p(1 / a)
    name == "LFMO-IG" && return (sqrt(2X + a^2) - a) / (sqrt(2 + a^2) - a)
    name == "LFMO-S" && return X^a
    name == "BFR" && return (1 - (1 - a)^X) / a
    name == "LSh" && return (1 - a) * X + a
    error("unknown model $(name); expected one of $(MODELS)")
end

"""
    psi(model, x) -> BigFloat

Laplace exponent Ψ(x) of the calibrated model, for integer x ≥ 0.
"""
psi(m::LFMOModel, x::Integer) = setprecision(() -> _psi_big(m.name, m.param, x), BigFloat, BIG_PRECISION)

"""
    lfmo_model(name, beta::AbstractString) -> LFMOModel

Calibrate the parameter of model `name` so that Ψ(1) = 1 and 2 - Ψ(2) = β.
"""
function lfmo_model(name::AbstractString, beta::AbstractString)
    name in MODELS || error("unknown model $(name); expected one of $(MODELS)")
    setprecision(BigFloat, BIG_PRECISION) do
        b = parse(BigFloat, beta)
        if name == "Exp"
            b == 0 || error("model Exp has β = 0, got β = $(beta)")
            return LFMOModel(name, beta, BigFloat(0))
        end
        0 < b < 1 || error("β must be in (0,1), got $(beta)")
        a = if name == "LFMO-C"
            b / (2 - b)
        elseif name in ("BFR", "LSh")
            b
        elseif name == "LFMO-S"
            log2(2 - b)
        else   # LFMO-G, LFMO-IG: β(η) = 2 - Ψ(2) is decreasing in η; bisection on log η
            g = eta -> 2 - _psi_big(name, eta, 2) - b
            lo, hi = BigFloat("1e-40"), BigFloat("1e40")
            g(lo) > 0 || error("β = $(beta) is not attainable by $(name) (sup β = $(Float64(2 - _psi_big(name, lo, 2))))")
            g(hi) < 0 || error("β = $(beta) too small for $(name)")
            for _ in 1:1000
                mid = sqrt(lo * hi)
                g(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            sqrt(lo * hi)
        end
        m = LFMOModel(name, beta, a)
        abs(_psi_big(name, a, 1) - 1) < BigFloat(10)^-100 || error("calibration failed: Ψ(1) ≠ 1 for $(name)")
        abs(2 - _psi_big(name, a, 2) - b) < BigFloat(10)^-100 || error("calibration failed: β ≠ $(beta) for $(name)")
        return m
    end
end

"""
    transition_rates(model, n) -> (R, exits)

Survivor-count chain of the LFMO model with `n` components, indexed by the
number `c` of failed components (array index c+1):

- `R[c+1, c2+1]` (c2 > c) is the rate at which the count jumps from c to c2,
  i.e. the rate at which exactly f = c2 - c of the w = n - c working
  components fail together, `binomial(w, f) λ^(w)_f` with `λ^(w)_f` as in
  eq. (def:lambda) of the paper;
- `exits[c+1] = Ψ(n - c)`, the total jump rate from c.

Computed in BigFloat and returned as Float64. Throws if a rate is negative or
the rates out of a state do not sum to Ψ(w) (both would indicate a Ψ that is
not a Bernstein function, or a loss of precision).
"""
function transition_rates(m::LFMOModel, n::Integer)
    setprecision(BigFloat, BIG_PRECISION) do
        P = [_psi_big(m.name, m.param, x) for x in 0:n]      # P[x+1] = Ψ(x)
        R = zeros(Float64, n + 1, n + 1)
        exits = zeros(Float64, n + 1)
        tol = BigFloat(10)^-80
        for c in 0:(n - 1)
            w = n - c
            total = BigFloat(0)
            for f in 1:w
                lam = BigFloat(0)
                for l in (w - f):w
                    term = binomial(big(f), w - l) * P[l + 1]
                    lam += isodd(l - (w - f) + 1) ? -term : term
                end
                rate = binomial(big(w), f) * lam
                rate > -tol * P[w + 1] || error("$(m.name) β=$(m.beta): negative rate $(Float64(rate)) (w=$(w), f=$(f))")
                rate = max(rate, BigFloat(0))
                total += rate
                R[c + 1, c + f + 1] = Float64(rate)
            end
            abs(total - P[w + 1]) <= tol * P[w + 1] || error("$(m.name) β=$(m.beta): rates out of w=$(w) do not sum to Ψ(w)")
            exits[c + 1] = Float64(P[w + 1])
        end
        return R, exits
    end
end

"""
    mean_batch(model, w) -> Float64

Expected number of components failing together at the next failure epoch when
`w` components work: w Ψ(1) / Ψ(w) = 1/κ(w).
"""
mean_batch(m::LFMOModel, w::Integer) = Float64(w * psi(m, 1) / psi(m, w))
