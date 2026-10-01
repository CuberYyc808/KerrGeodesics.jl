# Radial parts of t, φ and τ along an analytic r(λ).
#
#     t_r = ∫ (r² + a²) P/Δ dλ,   φ_r = ∫ (aP/Δ − aE) dλ = ∫ a(2Er − aLz)/Δ dλ,   τ_r = ∫ r² dλ
#
# The Mino-time domain is cut into segments, each integrated in the variable in which its
# integrand is smooth, and stored as a Chebyshev primitive (Spectral.jl):
#
#   :plain      λ      the rates themselves (no horizon, no infinity, no repeated root)
#   :horizon    λ      next to a horizon r_h (r₊ or r₋) the rates have a 1/(r − r_h) pole.
#                      With (P − s√R)/Δ = K/(P + s√R), K = r² + (Lz − aE)² + Q, s = sign P(r_h),
#                          (r²+a²)P/Δ = (r²+a²)K/(P+s√R) + sσ dr*/dλ,
#                          aP/Δ − aE  = aK/(P+s√R) − aE + sσ dφ_H/dλ,
#                      (σ = sign dr/dλ): the first terms are regular and are integrated,
#                      r* and φ_H (log|r − r±|) are added in closed form. v = t + r* (or
#                      u = t − r*) is then formed without cancellation. A segment may
#                      contain the horizon itself: t, φ are then continued through it as the
#                      principal value (t ∓ r* stays analytic across the crossing).
#   :infinity  q       an end at r → ∞ (finite λ), q = √(r_s/r) ∈ [0, 1] (r_s the split
#                      radius), q = 0 at infinity itself. With R = r⁴ W(q²),
#                      W(u) = c₄ + βu + γu² + δu³ + εu⁴ (β = 2/r_s, γ = c₂/r_s², …),
#                      dλ = −2q dq/(σ r_s √W): each rate becomes F(q)/√W with
#                      F = F₋₃ q⁻³ + F₋₁ q⁻¹ + F_reg (t: F₋₃ = −2E r_s/σ, F₋₁ = −4E/σ;
#                      τ: F₋₃ = −2r_s/σ; φ: regular). Writing W = S² + D, S² = c₄ + βq²,
#                          F/√W = (F₋₃ q⁻³ + F₋₁ q⁻¹)/S − (γ/2) F₋₃ q/S³ + rest(q)
#                      the first terms integrate in closed form (the r and log r terms
#                      for E > 1, the r^{3/2}, r^{1/2} terms for E = 1, and the passage
#                      between them for E → 1⁺ in one expression) and rest is bounded
#                      for every c₄ ≥ 0; only rest is fitted. Valid up to r = ∞.
#   :asymptote r       an end approaching a double root r_d as λ → ±∞:
#                          t_r = T_r(r_d) λ + σ ∫ (T_r − T_r(r_d))/√R dr,
#                      with √R = |r − r_d| √R₂ and the difference quotient taken from
#                      deflated polynomials, so r saturating at r_d costs nothing.
#              s       the same at a triple root: R₂ = (r − r_d) R₃ still vanishes at r_d,
#                      and the integrand keeps |r − r_d|^(−1/2); with r = r_d ± s² it is
#                      2σ Q(r)/√(±R₃) ds, smooth up to s = 0.
#
# Stable librations are periodic in λ and use one :plain period; a radial oscillation
# through the horizons (a plunge continued past r₊, r₋) is periodic over its legs'
# segments.

const _RADIAL_ENGINE_TOL = 1.0e-14

# ---- polynomials (ascending coefficients) ------------------------------------------------
@inline _horner(c, x) = (s = 0.0; for k in length(c):-1:1; s = muladd(s, x, c[k]); end; s)

function _deflate(c::Vector{Float64}, root)          # c(x) = (x − root) q(x) + rem
    n = length(c)
    n <= 1 && return Float64[]
    q = zeros(n - 1)
    acc = c[n]
    for k in n-1:-1:1
        q[k] = acc
        acc = c[k] + acc * root
    end
    return q
end

_polyaxpy(α, p, β, q) = [α * get(p, k, 0.0) + β * get(q, k, 0.0) for k in 1:max(length(p), length(q))]

# ---- rates --------------------------------------------------------------------------------
struct _RadialConstants{P}
    a::Float64; E::Float64; L::Float64; Q::Float64
    rplus::Float64
    R::Vector{Float64}          # radial potential coefficients (the infinity tail, the asymptote)
    potential::P                # R(r) for √R: the caller's form (product over the roots for the
                                # members, so R keeps its digits next to the horizon and next to
                                # a repeated root; the coefficient form for the frozen interfaces)
end

_rc(a, E, L, Q, potential) = _RadialConstants(a, E, L, Q, kerr_horizons(a).rplus,
    collect(Float64, kerr_radial_coefficients(a, E, L, Q)), potential)

@inline _rc_sqrtR(c::_RadialConstants, r) = sqrt(max(c.potential(r), 0.0))

# the coefficient form of R by Horner's rule (the APEX and finite-window interfaces)
function _coefficient_potential(a, E, L, Q)
    coefficients = collect(Float64, kerr_radial_coefficients(a, E, L, Q))
    return r -> _horner(coefficients, r)
end

@inline function _plain_rates(c::_RadialConstants, r)
    a, E, L = c.a, c.E, c.L
    Δ = kerr_delta(a, r)
    P = kerr_radial_momentum(a, E, L, r)
    return ((r^2 + a^2) * P / Δ, a * (2E * r - a * L) / Δ, r^2)
end

@inline function _horizon_rates(c::_RadialConstants, r, hs)
    a, E, L, Q = c.a, c.E, c.L, c.Q
    P = kerr_radial_momentum(a, E, L, r)
    K = r^2 + (L - a * E)^2 + Q
    D = P + hs * _rc_sqrtR(c, r)
    return ((r^2 + a^2) * K / D, a * K / D - a * E, r^2)
end

# The rounding of the radius moves the rates by eps·r·|∂rate/∂r|: the floor below which no
# fit of them can go. It matters for a leg that ends at a horizon that is nearly a root of R
# (D = P + s√R small; a plunge from a root a sliver above the horizon lives entirely there) and
# is negligible elsewhere. The r-derivatives in closed form, R' from the coefficients.
function _rate_rounding(kind, c::_RadialConstants, r, hs)
    a, E, L, Q = c.a, c.E, c.L, c.Q
    P = kerr_radial_momentum(a, E, L, r)
    if kind === :horizon
        K = r^2 + (L - a * E)^2 + Q
        rootR = _rc_sqrtR(c, r)
        D = P + hs * rootR
        Rprime = kerr_radial_derivatives(a, E, L, Q, r).R1
        Dprime = 2E * r + (iszero(rootR) ? 0.0 : hs * Rprime / (2rootR))
        dt = (2r * K + 2r * (r^2 + a^2)) / D - (r^2 + a^2) * K * Dprime / D^2
        dphi = 2r * a / D - a * K * Dprime / D^2
        return eps(Float64) * abs(r) .* (abs(dt), abs(dphi), 2abs(r))
    end
    Δ = kerr_delta(a, r); Δprime = 2 * (r - 1)
    dt = (2r * P + (r^2 + a^2) * 2E * r) / Δ - (r^2 + a^2) * P * Δprime / Δ^2
    dphi = 2E * a / Δ - a * (2E * r - a * L) * Δprime / Δ^2
    return eps(Float64) * abs(r) .* (abs(dt), abs(dphi), 2abs(r))
end

# ---- the end at infinity (header, :infinity) ------------------------------------------------
struct _InfinityTail
    rs::Float64                                  # split radius (q = 1)
    W::NTuple{5,Float64}                         # c₄, β, γ, δ, ε of W(u)
    F3::NTuple{3,Float64}; F1::NTuple{3,Float64} # coefficients of q⁻³, q⁻¹ (t, φ, τ)
    σ::Float64
end
const _NO_TAIL = _InfinityTail(NaN, (0.0, 0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0)

function _InfinityTail(c::_RadialConstants, rs, σ)
    c0, c1, c2, c3, c4 = c.R
    E = c.E
    return _InfinityTail(rs, (c4, c3 / rs, c2 / rs^2, c1 / rs^3, c0 / rs^4),
        (-2E * rs / σ, 0.0, -2rs / σ), (-4E / σ, 0.0, 0.0), σ)
end

# artanh(x) given x and 1 − x (1 − x carries the digits when x → 1)
_atanh_stable(x, omx) = x < 0.5 ? atanh(x) : 0.5 * log((2 - omx) / omx)
# (artanh x − x)/x³ = 1/3 + x²/5 + x⁴/7 + …
function _atanh_g1(x, omx)
    x < 0.1 || return (_atanh_stable(x, omx) - x) / x^3
    s = 0.0; p = 1.0
    for k in 1:12
        s += p / (2k + 1); p *= x^2
    end
    return s
end

# closed-form part of the tail primitives at q > 0 (zero at no particular point: only
# differences are used); diverges as q → 0 like the coordinates themselves
function _tail_principal(tl::_InfinityTail, q)
    c4, β, γ = tl.W[1], tl.W[2], tl.W[3]
    S = sqrt(c4 + β * q^2)
    κ = sqrt(c4)
    x = κ / S
    omx = β * q^2 / (S * (S + κ))                # 1 − x without cancellation
    I1 = -1 / (2q^2 * S) + β * _atanh_g1(x, omx) / (2S^3)   # ∫ dq/(q³ S)
    I2 = -(iszero(κ) ? 1.0 : _atanh_stable(x, omx) / x) / S # ∫ dq/(q S)
    I3 = -1 / (β * S)                                         # ∫ q dq/S³
    return ntuple(k -> tl.F3[k] * (I1 - γ / 2 * I3) + tl.F1[k] * I2, 3)
end

# the bounded rest of the tail rates (header), at q ∈ [0, 1]; q = 0 is r = ∞, where the
# expressions below are 0/0 for c₄ = 0 and take their limits
function _tail_rest(c::_RadialConstants, tl::_InfinityTail, q)
    c4, β, γ, δ, ε = tl.W
    a, E, L = c.a, c.E, c.L
    x = q^2 / tl.rs                              # 1/r
    Δx = 1 - 2x + a^2 * x^2                      # Δ/r²
    # regular parts of the rates: t − E r² − 2E r, φ (both finite at r = ∞)
    Treg = ((E * a^2 - a * L + 4E) - 2E * a^2 * x + a^2 * (E * a^2 - a * L) * x^2) / Δx
    Φ = a * (2E - a * L * x) * x / Δx
    Dt = γ + δ * q^2 + ε * q^4                   # D/q⁴
    if iszero(q)
        iszero(c4) || return (0.0, 0.0, 0.0)
        s = sqrt(β)                              # S/q = √W/q at q = 0
        k3 = 3γ^2 / (8s^5) - δ / (2s^3)
        k1 = γ / (2s^3)
        w = -2 / (tl.σ * tl.rs * s)
        return (w * Treg + tl.F3[1] * k3 - tl.F1[1] * k1, w * Φ, tl.F3[3] * k3)
    end
    S = sqrt(c4 + β * q^2)
    rW = sqrt(S^2 + Dt * q^4)
    D = Dt * q^4
    den = rW * S * (S + rW)
    k3 = γ * q * D * (2S + rW) / (2S^3 * rW * (S + rW)^2) - (δ + ε * q^2) * q^3 / den
    k1 = Dt * q^3 / den
    w = -2q / (tl.σ * tl.rs * rW)
    return (w * Treg + tl.F3[1] * k3 - tl.F1[1] * k1, w * Φ, tl.F3[3] * k3)
end

_tail_q(tl::_InfinityTail, r) = sqrt(tl.rs / r)

# ---- segments -----------------------------------------------------------------------------
mutable struct _RadialSegment
    kind::Symbol
    lo::Float64; hi::Float64                 # λ range
    σ::Float64                               # sign of dr/dλ
    anchor_λ::Float64                        # λ where the segment's value equals `base`
    prim::ChebPieces
    prim_anchor::NTuple{3,Float64}
    base::NTuple{3,Float64}
    rd::Float64                              # :asymptote only
    multiplicity::Int                        # :asymptote only: 2, or 3 (variable s)
    rates_d::NTuple{3,Float64}
    tail::_InfinityTail                      # :infinity only
    rstar_anchor::Float64; azimuth_anchor::Float64  # :horizon only
    hs::Float64                                     # :horizon only: sign P(r_h)
end

_needs_r(s::_RadialSegment) = s.kind !== :plain

function _segment_value(s::_RadialSegment, c::_RadialConstants, λ, r)
    kind = s.kind
    if kind === :plain
        return s.base .+ (_eval3(s.prim, float(λ)) .- s.prim_anchor)
    elseif kind === :horizon
        F = _eval3(s.prim, float(λ)) .- s.prim_anchor
        w = s.σ * s.hs
        return (s.base[1] + F[1] + w * (_rstar_all(c.a, r) - s.rstar_anchor),
                s.base[2] + F[2] + w * (_horizon_azimuth(c.a, r) - s.azimuth_anchor),
                s.base[3] + F[3])
    elseif kind === :asymptote
        return s.base .+ (_eval3(s.prim, _asymptote_x(s, r)) .- s.prim_anchor) .+
            s.rates_d .* (λ - s.anchor_λ)
    end
    return s.base .+ (_tail_value(s, r) .- s.prim_anchor)          # :infinity
end

# the fit variable of an :asymptote segment at radius r (r, or s = √|r − r_d|)
function _asymptote_x(s::_RadialSegment, r)
    br = s.prim.breaks
    x = s.multiplicity == 3 ? sqrt(abs(r - s.rd)) : r
    return clamp(x, min(br[1], br[end]), max(br[1], br[end]))
end

function _tail_value(s::_RadialSegment, r)
    q = min(_tail_q(s.tail, r), 1.0)
    return _tail_principal(s.tail, q) .+ _eval3(s.prim, q)
end

# the horizon-regular parts: t − sσ r*, φ − sσ φ_H (finite at the horizon this segment meets)
function _segment_regular(s::_RadialSegment, λ)
    F = _eval3(s.prim, float(λ)) .- s.prim_anchor
    w = s.σ * s.hs
    return (s.base[1] + F[1] - w * s.rstar_anchor, s.base[2] + F[2] - w * s.azimuth_anchor)
end

function _build_segment(kind, c::_RadialConstants, r_of, lo, hi, σ; rd=NaN, multiplicity=2,
        rh=c.rplus)
    a = c.a; E = c.E
    floors = (0.0, abs(a) * (1 + abs(E)), 0.0)
    rates_d = (0.0, 0.0, 0.0); tail = _NO_TAIL; hs = 1.0
    if kind === :plain
        f = λ -> _plain_rates(c, r_of(λ))
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=_RADIAL_ENGINE_TOL, absfloor=floors,
            abserr=λ -> _rate_rounding(:plain, c, r_of(λ), hs)))
    elseif kind === :horizon
        hs = kerr_radial_momentum(a, E, c.L, rh) >= 0 ? 1.0 : -1.0
        f = λ -> _horizon_rates(c, r_of(λ), hs)
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=_RADIAL_ENGINE_TOL, absfloor=floors,
            abserr=λ -> _rate_rounding(:horizon, c, r_of(λ), hs)))
    elseif kind === :infinity
        # variable q = √(r_s/r) from the finite end (q = 1) to infinity (q = 0)
        tail = _InfinityTail(c, r_of(σ > 0 ? lo : hi), σ)
        prim = chebintegrate(chebfit(q -> _tail_rest(c, tail, q), 0.0, 1.0; ncomp=3,
            tol=_RADIAL_ENGINE_TOL, absfloor=floors))
    elseif kind === :asymptote
        aa, EE, LL = c.a, c.E, c.L
        num_t = [EE * aa^4 - aa * LL * aa^2, 0.0, 2EE * aa^2 - aa * LL, 0.0, EE]
        num_φ = aa .* [2EE * 0.0 - aa * LL, 2EE]           # a(2E r − aL)
        Δc = [aa^2, -2.0, 1.0]
        Δd = _horner(Δc, rd)
        Qt = _deflate(_polyaxpy(Δd, num_t, -_horner(num_t, rd), Δc), rd)
        Qφ = _deflate(_polyaxpy(Δd, num_φ, -_horner(num_φ, rd), Δc), rd)
        R2 = _deflate(_deflate(c.R, rd), rd)
        rates_d = (_horner(num_t, rd) / Δd, _horner(num_φ, rd) / Δd, rd^2)
        rfar = r_of(isfinite(lo) ? lo : hi)                # the finite end
        sgn = rfar >= rd ? 1.0 : -1.0                      # side of r_d the orbit is on
        if multiplicity == 3
            R3 = _deflate(R2, rd)
            g = function (x)
                r = rd + sgn * x^2
                w = 2σ / sqrt(sgn * _horner(R3, r))
                Δr = _horner(Δc, r)
                return (w * _horner(Qt, r) / (Δr * Δd), w * _horner(Qφ, r) / (Δr * Δd),
                    w * (r + rd))
            end
            prim = chebintegrate(chebfit(g, 0.0, sqrt(abs(rfar - rd)); ncomp=3,
                tol=_RADIAL_ENGINE_TOL, absfloor=floors))
        else
            f = function (r)
                w = σ * sgn / sqrt(max(_horner(R2, r), floatmin()))
                Δr = _horner(Δc, r)
                return (w * _horner(Qt, r) / (Δr * Δd), w * _horner(Qφ, r) / (Δr * Δd),
                    w * (r + rd))
            end
            prim = chebintegrate(chebfit(f, min(rd, rfar), max(rd, rfar); ncomp=3,
                tol=_RADIAL_ENGINE_TOL, absfloor=floors))
        end
    else
        error("unknown radial segment kind $kind")
    end
    return _RadialSegment(kind, lo, hi, σ, NaN, prim, (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0), rd, multiplicity, rates_d, tail, NaN, NaN, hs)
end

function _anchor!(s::_RadialSegment, c, r_of, λ, base)
    r = _needs_r(s) ? r_of(λ) : NaN
    s.anchor_λ = λ
    s.prim_anchor = s.kind === :infinity ? _tail_value(s, r) :
        _eval3(s.prim, s.kind === :asymptote ? _asymptote_x(s, r) : λ)
    s.base = base
    if s.kind === :horizon
        s.rstar_anchor = _rstar_all(c.a, r)
        s.azimuth_anchor = _horizon_azimuth(c.a, r)
    end
    return s
end

# λ on a monotone leg where r(λ) = target: bisection between an inner point and the far
# end (a finite endpoint is never evaluated; an infinite one is bracketed by doubling).
function _leg_lambda(r_of, λin, λout, target)
    rin = r_of(λin)
    dir = sign(λout - λin)
    if isfinite(λout)
        lo, hi = λin, λout
    else
        step = 1.0
        lo = λin; hi = λin + dir * step
        while (r_of(hi) - target) * (rin - target) > 0
            lo = hi; step *= 2; hi = λin + dir * step
            abs(step) > 1e6 && error("radial engine: cannot bracket r = $target")
        end
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        if (r_of(mid) - target) * (rin - target) > 0
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

"""
    RadialEngine

Radial primitives along r(λ) on a Mino-time domain, zero at `λ_ref`. Build with
`_radial_engine`; evaluate with `_radial_eval(e, λ)` → (t_r, φ_r, τ_r) and
`_radial_regular(e, λ, σ)` → (t_r − σ r*, φ_r − σ φ_H), the combination that stays finite
at the horizon the orbit meets (σ = −1: v, ψ; σ = +1: u, χ).
"""
struct RadialEngine{F,P}
    c::_RadialConstants{P}
    r_of::F
    segments::Vector{_RadialSegment}
    periodic::Bool
    period::Float64
    totals::NTuple{3,Float64}
    shift::NTuple{3,Float64}
end

"""
    _radial_engine(a, E, Lz, Q, r_of; potential, domain, ends, turn=nothing, σ, rd=NaN,
                   multiplicity=2, λ_ref=0.0)

`potential(r)` is R(r) as the caller forms it (see `_RadialConstants`).

`ends = (kind_lo, kind_hi)` with kinds `:turning` (any finite regular point), `:horizon`, `:infinity`,
`:asymptote`; `turn` is an interior turning point (or `nothing`) and `σ` the sign of
dr/dλ just above `turn` (or on the whole domain); `rd` the repeated root of an asymptote and
`multiplicity` its order (2 or 3).
For a periodic libration pass `period` and `domain = (λ_peri, λ_peri + period)`; with a
`turn` as well (ends are then turning points) the period is cut into its two legs, which
may cross the horizons.
"""
function _radial_engine(a, E, L, Q, r_of; potential, domain, ends=(:turning, :turning),
        turn=nothing, σ=1.0, rd=NaN, multiplicity=2, λ_ref=0.0, period=nothing)
    t = _radial_engine_tables(float(a), float(E), float(L), float(Q), _ErasedFunction(r_of),
        _ErasedFunction(potential), (float(domain[1]), float(domain[2])), (ends[1], ends[2]),
        turn === nothing ? nothing : float(turn), float(σ), float(rd), Int(multiplicity),
        float(λ_ref), period === nothing ? nothing : float(period))
    return RadialEngine(_rc(a, E, L, Q, potential), r_of, t.segments, t.periodic, t.period,
        t.totals, t.shift)
end

# A radius or potential function behind a field of abstract type. The engine's tables are built
# through it, so the table construction (segment search, Chebyshev fits, anchoring) is compiled
# once for every radial model instead of once per model, at the price of a dynamic call per
# sample; the engine evaluates the tables with the concrete functions.
struct _ErasedFunction
    f::Any
end
(e::_ErasedFunction)(x) = convert(Float64, e.f(x))::Float64

function _radial_engine_tables(a, E, L, Q, r_of::_ErasedFunction, potential::_ErasedFunction,
        domain, ends, turn, σ, rd, multiplicity, λ_ref, period)
    c = _rc(a, E, L, Q, potential)
    if period !== nothing && turn !== nothing
        e = _radial_engine_tables(a, E, L, Q, r_of, potential, domain, (:turning, :turning), turn,
            σ, NaN, 2, domain[1], nothing)
        lo = float(domain[1])
        totals = _radial_eval_raw(e, lo + period)
        e = RadialEngine(c, r_of, e.segments, true, float(period), totals, (0.0, 0.0, 0.0))
        v = _radial_eval_raw(e, float(λ_ref))
        return RadialEngine(c, r_of, e.segments, true, float(period), totals, v)
    elseif period !== nothing
        lo = domain[1]
        s = _build_segment(:plain, c, r_of, lo, lo + period, 1.0)
        _anchor!(s, c, r_of, lo, (0.0, 0.0, 0.0))
        totals = ntuple(k -> s.prim(lo + period, k) - s.prim_anchor[k], 3)
        e = RadialEngine(c, r_of, [s], true, float(period), totals, (0.0, 0.0, 0.0))
        v = _radial_eval_raw(e, float(λ_ref))
        return RadialEngine(c, r_of, [s], true, float(period), totals, v)
    end
    λa, λb = float.(domain)
    legs = turn === nothing ? [(λa, λb, ends[1], ends[2], float(σ))] :
        [(λa, float(turn), ends[1], :turning, -float(σ)), (float(turn), λb, :turning, ends[2], float(σ))]
    segs = _RadialSegment[]
    for (lo, hi, klo, khi, sg) in legs
        append!(segs, _leg_segments(c, r_of, lo, hi, klo, khi, sg, rd, multiplicity, λ_ref))
    end
    sort!(segs; by=s -> s.lo)
    # anchor: the segment holding λ_ref is zero there; neighbours continue from it
    i0 = findfirst(s -> s.lo <= λ_ref <= s.hi, segs)
    i0 === nothing && error("radial engine: reference λ = $λ_ref outside the domain.")
    _anchor!(segs[i0], c, r_of, λ_ref, (0.0, 0.0, 0.0))
    for j in i0+1:length(segs)
        b = segs[j].lo
        _anchor!(segs[j], c, r_of, b, _segment_value(segs[j-1], c, b, r_of(b)))
    end
    for j in i0-1:-1:1
        b = segs[j].hi
        _anchor!(segs[j], c, r_of, b, _segment_value(segs[j+1], c, b, r_of(b)))
    end
    return RadialEngine(c, r_of, segs, false, 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
end

_end_radius(c, r_of, λ, kind) = kind === :horizon ? c.rplus : r_of(λ)

function _leg_segments(c, r_of, lo, hi, klo, khi, σ, rd, multiplicity, λ_ref)
    segs = _RadialSegment[]
    open_end(k) = k === :infinity || k === :asymptote      # ends that get their own segment
    # finite inner reference radius of the leg
    rin_lo = open_end(klo) ? NaN : _end_radius(c, r_of, lo, klo)
    rin_hi = open_end(khi) ? NaN : _end_radius(c, r_of, hi, khi)
    # an interior point of the leg where r is finite and regular
    nudge = 1e-9 * (isfinite(hi - lo) ? max(1.0, abs(hi - lo)) : 1.0)
    inner = !open_end(klo) ? lo + nudge : !open_end(khi) ? hi - nudge : float(λ_ref)
    rref = isfinite(rin_lo) ? rin_lo : isfinite(rin_hi) ? rin_hi : r_of(inner)
    function split_for(k, λend)
        if k === :infinity
            target = max(2rref, 20.0)
        else
            target = rd + 0.5 * (rref - rd)
        end
        return _leg_lambda(r_of, inner, λend, target)
    end
    s_lo = open_end(klo) ? split_for(klo, lo) : lo
    s_hi = open_end(khi) ? split_for(khi, hi) : hi
    open_end(klo) && push!(segs, _build_segment(klo, c, r_of, lo, s_lo, σ; rd=rd,
        multiplicity=multiplicity))
    # the finite middle: horizon kernel next to a horizon, plain rates elsewhere
    hor_lo = klo === :horizon; hor_hi = khi === :horizon
    crossed = open_end(klo) || open_end(khi) || hor_lo || hor_hi ? Float64[] :
        _crossed_horizons(c, r_of(s_lo), r_of(s_hi))
    if !isempty(crossed)
        # turning/regular ends with horizons in between: one horizon segment around each
        # crossing (continued through it), plain rates next to the ends; cuts at the
        # Mino-time midpoints between consecutive crossing points
        λh = [_leg_lambda(r_of, s_lo, s_hi, h) for h in crossed]
        pts = [s_lo; λh; s_hi]
        cuts = [(pts[k] + pts[k+1]) / 2 for k in 1:length(pts)-1]
        push!(segs, _build_segment(:plain, c, r_of, s_lo, cuts[1], σ))
        for (k, h) in enumerate(crossed)
            push!(segs, _build_segment(:horizon, c, r_of, cuts[k], cuts[k+1], σ; rh=h))
        end
        push!(segs, _build_segment(:plain, c, r_of, cuts[end], s_hi, σ))
    elseif hor_lo || hor_hi
        turning_other = (hor_lo ? khi : klo) === :turning
        if turning_other
            # A finite apastron can be arbitrarily far away as E approaches one.
            # Geometric radial pieces keep each fit local and avoid subtracting
            # a far-apastron primitive to obtain a small near-horizon increment.
            lambda_h=hor_lo ? s_lo : s_hi
            lambda_t=hor_lo ? s_hi : s_lo
            rt=r_of(lambda_t)
            rh=c.rplus
            count=max(2,ceil(Int,log2(rt/rh)))
            cuts=Float64[lambda_h]
            for j in 1:count-1
                target=rh*exp(log(rt/rh)*j/count)
                push!(cuts,_leg_lambda(r_of,lambda_h,lambda_t,target))
            end
            push!(cuts,lambda_t)
            for j in 1:count
                left,right=minmax(cuts[j],cuts[j+1])
                push!(segs,_build_segment(j==1 ? :horizon : :plain,c,r_of,left,right,σ))
            end
        else
            push!(segs, _build_segment(:horizon, c, r_of, s_lo, s_hi, σ))
        end
    else
        push!(segs, _build_segment(:plain, c, r_of, s_lo, s_hi, σ))
    end
    open_end(khi) && push!(segs, _build_segment(khi, c, r_of, s_hi, hi, σ; rd=rd,
        multiplicity=multiplicity))
    return segs
end

# horizons strictly between two radii, in the order met going from `ra` to `rb` (r₋ only for
# a ≠ 0, where it is a pole of the rates)
function _crossed_horizons(c::_RadialConstants, ra, rb)
    lo, hi = minmax(ra, rb)
    h = kerr_horizons(c.a)
    hs = Float64[x for x in (h.rplus, h.rminus) if lo < x < hi && !(x == h.rminus && iszero(c.a))]
    unique!(hs)
    return ra > rb ? sort!(hs; rev=true) : sort!(hs)
end

@inline function _find_segment(e::RadialEngine, λ)
    segs = e.segments
    @inbounds for i in eachindex(segs)
        λ <= segs[i].hi && return segs[i]
    end
    return segs[end]
end

function _radial_eval_raw(e::RadialEngine, λ)
    if e.periodic
        s = e.segments[1]
        n = floor((λ - s.lo) / e.period)
        if length(e.segments) == 1
            # a libration from periapsis: the rates are symmetric about apoapsis, so the second
            # half-period is integrated back from the next periapsis. (Adding a full period
            # to −totals would leave only eps·|totals| of absolute precision next to λ_peri.)
            x = λ - n * e.period - s.lo
            2x <= e.period && return (_eval3(s.prim, s.lo + x) .- s.prim_anchor) .+ n .* e.totals
            return (n + 1) .* e.totals .- (_eval3(s.prim, s.lo + (e.period - x)) .- s.prim_anchor)
        end
        μ = λ - n * e.period
        s = _find_segment(e, μ)
        return _segment_value(s, e.c, μ, _needs_r(s) ? e.r_of(μ) : NaN) .+ n .* e.totals
    end
    s = _find_segment(e, λ)
    return _segment_value(s, e.c, λ, _needs_r(s) ? e.r_of(λ) : NaN)
end
function _radial_eval_raw(e::RadialEngine, λ, r)
    if e.periodic
        length(e.segments) == 1 && return _radial_eval_raw(e, λ)
        n = floor((λ - e.segments[1].lo) / e.period)
        μ = λ - n * e.period
        return _segment_value(_find_segment(e, μ), e.c, μ, r) .+ n .* e.totals
    end
    return _segment_value(_find_segment(e, λ), e.c, λ, r)
end

"""(t_r, φ_r, τ_r) at λ (zero at the engine's reference λ); pass r if already known."""
_radial_eval(e::RadialEngine, λ) = _radial_eval_raw(e, λ) .- e.shift
_radial_eval(e::RadialEngine, λ, r) = _radial_eval_raw(e, λ, r) .- e.shift

function _radial_proper_interval(e::RadialEngine, left, right)
    left==right && return 0.0
    right<left && return -_radial_proper_interval(e,right,left)
    value=0.0
    for s in e.segments
        lo=max(left,s.lo); hi=min(right,s.hi)
        hi<=lo && continue
        if s.kind===:plain || s.kind===:horizon
            value+=_cheb_increment(s.prim,lo,hi,3)
        elseif s.kind===:asymptote
            xlo=_asymptote_x(s,e.r_of(lo)); xhi=_asymptote_x(s,e.r_of(hi))
            value+=_cheb_increment(s.prim,xlo,xhi,3)+s.rates_d[3]*(hi-lo)
        else
            qlo=min(_tail_q(s.tail,e.r_of(lo)),1.0)
            qhi=min(_tail_q(s.tail,e.r_of(hi)),1.0)
            value+=_tail_principal(s.tail,qhi)[3]-_tail_principal(s.tail,qlo)[3]+
                _cheb_increment(s.prim,qlo,qhi,3)
        end
    end
    return value
end

function _radial_proper_increment(e::RadialEngine, left, right)
    right<left && return -_radial_proper_increment(e,right,left)
    e.periodic || return _radial_proper_interval(e,left,right)
    lo=e.segments[1].lo
    nl=floor((left-lo)/e.period); nr=floor((right-lo)/e.period)
    nl==nr && return _radial_proper_interval(e,left-nl*e.period,right-nr*e.period)
    return (nr-nl-1)*e.totals[3]+
        _radial_proper_interval(e,left-nl*e.period,lo+e.period)+
        _radial_proper_interval(e,lo,right-nr*e.period)
end

"""
(t_r − σ r*, φ_r − σ φ_H): with σ = −1 the ingoing combination (v, ψ), with σ = +1 the
outgoing one (u, χ). Evaluated without cancellation on a horizon segment of that
orientation, so it stays finite at the horizon itself.
"""
function _radial_regular(e::RadialEngine, λ, σ_wanted)
    n = e.periodic && length(e.segments) > 1 ? floor((λ - e.segments[1].lo) / e.period) : 0.0
    μ = λ - n * e.period
    s = _find_segment(e, μ)
    if s.kind === :horizon && s.σ * s.hs == σ_wanted
        v = _segment_regular(s, μ)
        return (v[1] + n * e.totals[1] - e.shift[1], v[2] + n * e.totals[2] - e.shift[2])
    end
    r = e.r_of(λ)
    t, φ, _ = _radial_eval(e, λ, r)
    return (t - σ_wanted * kerr_rstar(e.c.a, r), φ - σ_wanted * _horizon_azimuth(e.c.a, r))
end

_spectral_summary(e::RadialEngine) = _spectral_summary(s.prim for s in e.segments)

"""Mean radial rates over one period (periodic engines)."""
_radial_mean_rates(e::RadialEngine) = e.totals ./ e.period

"""
    EngineCoordinates

The coordinates of one member: radial engine part + polar primitive (`polar_primitive(λ)`
gives the polar (t, φ, τ) increments). One object holds, once, the radius function, the polar
primitive, the potential, the engine's build parameters and the member's two lazy caches (the
radial engine, built on the first t/φ/τ call, and the origin of the regular chart), so that the
closures of a member capture one small object instead of trees of closures. Evaluate it with
`_coords_t`, `_coords_phi`, `_coords_tau`, `_coords_tphitau`, `_coords_v`, `_coords_psi`,
`_coords_radial`, `_coords_spectral`; `_regular_chart(c, σ, λ_ref)` gives further regular charts.
"""
struct EngineCoordinates{F,P,V}
    a::Float64; E::Float64; L::Float64; Q::Float64
    r_of::F
    polar_primitive::P
    potential::V
    domain::Tuple{Float64,Float64}
    ends::Tuple{Symbol,Symbol}
    turn::Union{Nothing,Float64}
    σ::Float64
    rd::Float64
    multiplicity::Int
    period::Union{Nothing,Float64}
    λ_bl::Float64                   # t = φ = 0 here
    pb::NTuple{3,Float64}           # polar primitive at λ_bl
    λ_tau::Float64                  # τ = 0 here
    pτ::Float64                     # polar τ primitive at λ_tau
    λ_regular::Union{Nothing,Float64}   # v = ψ = 0 here (nothing: no regular chart)
    σ_regular::Float64
    cache::Base.RefValue{Union{Nothing,RadialEngine{F,V}}}
    regular_origin::Base.RefValue{NTuple{2,Float64}}
end

"""
    _engine_coordinates(a, E, Lz, Q, r_of, polar_primitive; potential, domain, ends,
                        turn=nothing, σ=1.0, rd=NaN, multiplicity=2, λ_bl=0.0, λ_tau=0.0,
                        λ_regular=nothing, σ_regular=-1.0, period=nothing)

The `EngineCoordinates` of one member. `t`, `phi` vanish at `λ_bl`, `tau` at `λ_tau` (the
λ = 0 event of every member), and the horizon-regular pair `v = t − σ_regular r*`,
`psi = φ − σ_regular φ_H` at `λ_regular` (the horizon event, where they stay finite). The
engine is built on the first call.
"""
function _engine_coordinates(a, E, L, Q, r_of::F, polar_primitive::P; potential::V, domain,
        ends, turn=nothing, σ=1.0, rd=NaN, multiplicity=2, λ_bl=0.0, λ_tau=0.0,
        λ_regular=nothing, σ_regular=-1.0, period=nothing) where {F,P,V}
    return EngineCoordinates{F,P,V}(a, E, L, Q, r_of, polar_primitive, potential,
        (float(domain[1]), float(domain[2])), (ends[1], ends[2]),
        turn === nothing ? nothing : float(turn), σ, rd, multiplicity,
        period === nothing ? nothing : float(period), λ_bl, polar_primitive(λ_bl), λ_tau,
        polar_primitive(λ_tau)[3], λ_regular === nothing ? nothing : float(λ_regular),
        σ_regular, Ref{Union{Nothing,RadialEngine{F,V}}}(nothing), Ref((NaN, NaN)))
end

function _coords_engine(c::EngineCoordinates)
    e = c.cache[]
    e === nothing || return e
    e = _radial_engine(c.a, c.E, c.L, c.Q, c.r_of; potential=c.potential, domain=c.domain,
        ends=c.ends, turn=c.turn, σ=c.σ, rd=c.rd, multiplicity=c.multiplicity, λ_ref=c.λ_bl,
        period=c.period)
    c.cache[] = e
    return e
end

_coords_tphitau(c::EngineCoordinates, λ) =
    _radial_eval(_coords_engine(c), λ) .+ (c.polar_primitive(λ) .- c.pb)
_coords_t(c::EngineCoordinates, λ) = _coords_tphitau(c, λ)[1]
_coords_phi(c::EngineCoordinates, λ) = _coords_tphitau(c, λ)[2]
_coords_tau(c::EngineCoordinates, λ) =
    _radial_proper_increment(_coords_engine(c), c.λ_tau, λ) + c.polar_primitive(λ)[3] - c.pτ
_coords_radial(c::EngineCoordinates, λ) = _radial_eval(_coords_engine(c), λ)
_coords_spectral(c::EngineCoordinates) = _spectral_summary(_coords_engine(c))

# horizon-regular pair (t − σr r*, φ − σr φ_H) + polar part, zero at λr (origin cached in `origin`)
function _coords_regular(c::EngineCoordinates, σr, λr, origin, λ)
    e = _coords_engine(c)
    if isnan(origin[][1])
        v0 = _radial_regular(e, λr, σr)
        p0 = c.polar_primitive(λr)
        origin[] = (v0[1] + p0[1], v0[2] + p0[2])
    end
    v = _radial_regular(e, λ, σr)
    p = c.polar_primitive(λ)
    return (v[1] + p[1] - origin[][1], v[2] + p[2] - origin[][2])
end
function _coords_regular(c::EngineCoordinates, λ)
    c.λ_regular === nothing && error("No horizon-regular reference event for this member.")
    return _coords_regular(c, c.σ_regular, c.λ_regular, c.regular_origin, λ)
end
_coords_v(c::EngineCoordinates, λ) = _coords_regular(c, λ)[1]
_coords_psi(c::EngineCoordinates, λ) = _coords_regular(c, λ)[2]

"""A further regular chart (t − σr r*, φ − σr φ_H), zero at λr, with its own origin cache."""
function _regular_chart(c::EngineCoordinates, σr, λr)
    origin = Ref((NaN, NaN))
    return λ -> _coords_regular(c, σr, λr, origin, λ)
end

"""
    _radius_increments(coords, λ_of, σ; σ_regular=-1.0, regular=true)

Radial integrals between two radii of one monotone leg (dr/dλ = σ√R, `λ_of(r)` its Mino
time), ∫_{r1}^{r2} rate dr/√R, read off the engine of `_engine_coordinates`: the trajectory
fields `radial_mino/time/phi/proper_increment(r1, r2)` and, with `regular`, the
horizon-regular `radial_v/psi_increment` (t − σ_regular r*, φ − σ_regular φ_H, formed without
cancellation next to the horizon).
"""
function _radius_increments(coords, λ_of, σ; σ_regular=-1.0, regular=true)
    # the known radius is passed on: next to the horizon t, φ carry log|r − r₊|
    radial(k) = (r1, r2) -> σ * (_radial_eval(_coords_engine(coords), λ_of(r2), float(r2))[k] -
                                 _radial_eval(_coords_engine(coords), λ_of(r1), float(r1))[k])
    chart(k) = (r1, r2) -> σ * (_radial_regular(_coords_engine(coords), λ_of(r2), σ_regular)[k] -
                                _radial_regular(_coords_engine(coords), λ_of(r1), σ_regular)[k])
    increments = (radial_mino_increment=(r1, r2) -> σ * (λ_of(r2) - λ_of(r1)),
        radial_time_increment=radial(1), radial_phi_increment=radial(2),
        radial_proper_increment=radial(3))
    return regular ? merge(increments,
        (radial_v_increment=chart(1), radial_psi_increment=chart(2))) : increments
end
