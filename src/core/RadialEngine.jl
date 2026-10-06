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

# the relative accuracy of every radial table: 1e-14 in Float64, carried to T by `_tol`
_radial_engine_tol(::Type{T}) where {T} = _tol(T, 1.0e-14)

# ---- polynomials (ascending coefficients) ------------------------------------------------
@inline _horner(c, x) = (s = zero(promote_type(eltype(c), typeof(x))); for k in length(c):-1:1;
    s = muladd(s, x, c[k]); end; s)

function _deflate(c::Vector{T}, root) where {T}      # c(x) = (x − root) q(x) + rem
    n = length(c)
    n <= 1 && return T[]
    q = zeros(T, n - 1)
    acc = c[n]
    for k in n-1:-1:1
        q[k] = acc
        acc = c[k] + acc * root
    end
    return q
end

_polyaxpy(α, p, β, q) = [α * get(p, k, zero(eltype(p))) + β * get(q, k, zero(eltype(q)))
    for k in 1:max(length(p), length(q))]

# ---- rates --------------------------------------------------------------------------------
struct _RadialConstants{T,P}
    a::T; E::T; L::T; Q::T
    rplus::T
    R::Vector{T}                # radial potential coefficients (the infinity tail, the asymptote)
    potential::P                # R(r) for √R: the caller's form (product over the roots for the
                                # members, so R keeps its digits next to the horizon and next to
                                # a repeated root; the coefficient form for the frozen interfaces)
end

function _rc(a, E, L, Q, potential)
    T = _float_type(a, E, L, Q)
    return _RadialConstants{T,typeof(potential)}(a, E, L, Q, kerr_horizons(T(a)).rplus,
        collect(T, kerr_radial_coefficients(a, E, L, Q)), potential)
end

@inline _rc_sqrtR(c::_RadialConstants, r) = sqrt(max(c.potential(r), 0.0))

# Keep the small P near a horizon without rounding E(1+a^2)-aL first.
@noinline function _rc_momentum(a, energy, lz, r)
    a,E,L,x = _wide.(float.((a,energy,lz,r)))
    square = _wide_add(_wide_mul(x,x),_wide_mul(a,a))
    p = _wide_sub(_wide_mul(E,square),_wide_mul(a,L))
    return p[1]+p[2]
end
@inline _rc_momentum(c::_RadialConstants,r) = _rc_momentum(c.a,c.E,c.L,r)

@noinline function _rc_azimuth_numerator(a, energy, lz, r)
    a,E,L,x = _wide.(float.((a,energy,lz,r)))
    n = _wide_sub(_wide_mul(_wide(2.0),_wide_mul(E,x)),_wide_mul(a,L))
    return n[1]+n[2]
end
@inline _rc_azimuth_numerator(c::_RadialConstants,r) =
    _rc_azimuth_numerator(c.a,c.E,c.L,r)

_radial_state(r_of,lambda) = nothing
function _engine_rates(kind,c,r_of,lambda,hs)
    state = _radial_state(r_of,lambda)
    state === nothing || return _relative_rates(c,state,kind,hs)
    radius = r_of(lambda)
    return kind === :horizon ? _horizon_rates(c,radius,hs) : _plain_rates(c,radius)
end
function _engine_rstar(c,r_of,lambda,radius)
    state = _radial_state(r_of,lambda)
    return state === nothing ? _rstar_all(c.a,radius) : _relative_rstar(state)
end
function _engine_azimuth(c,r_of,lambda,radius)
    state = _radial_state(r_of,lambda)
    return state === nothing ? _horizon_azimuth(c.a,radius) : _relative_azimuth(c.a,state)
end
function _engine_rounding(kind,c::_RadialConstants{T},r_of,lambda,hs) where {T}
    _radial_state(r_of,lambda) === nothing || return (zero(T),zero(T),zero(T))
    return _rate_rounding(kind,c,r_of(lambda),hs)
end

# the coefficient form of R by Horner's rule (the APEX and finite-window interfaces)
function _coefficient_potential(a, E, L, Q)
    coefficients = collect(_float_type(a, E, L, Q), kerr_radial_coefficients(a, E, L, Q))
    return r -> _horner(coefficients, r)
end

@inline function _plain_rates(c::_RadialConstants, r)
    a, E, L = c.a, c.E, c.L
    Δ = kerr_delta(a, r)
    P = _rc_momentum(c,r)
    return ((r^2 + a^2) * P / Δ, a * _rc_azimuth_numerator(c,r) / Δ, r^2)
end

@inline function _horizon_rates(c::_RadialConstants, r, hs)
    a, E, L, Q = c.a, c.E, c.L, c.Q
    P = _rc_momentum(c,r)
    K = r^2 + (L - a * E)^2 + Q
    D = P + hs * _rc_sqrtR(c, r)
    return ((r^2 + a^2) * K / D, a * K / D - a * E, r^2)
end

# The rounding of the radius moves the rates by eps·r·|∂rate/∂r|: the floor below which no
# fit of them can go. It matters for a leg that ends at a horizon that is nearly a root of R
# (D = P + s√R small; a plunge from a root a sliver above the horizon lives entirely there) and
# is negligible elsewhere. The r-derivatives in closed form, R' from the coefficients.
function _rate_rounding(kind, c::_RadialConstants{T}, r, hs) where {T}
    a, E, L, Q = c.a, c.E, c.L, c.Q
    P = _rc_momentum(c,r)
    if kind === :horizon
        K = r^2 + (L - a * E)^2 + Q
        rootR = _rc_sqrtR(c, r)
        D = P + hs * rootR
        Rprime = kerr_radial_derivatives(a, E, L, Q, r).R1
        Dprime = 2E * r + (iszero(rootR) ? zero(T) : hs * Rprime / (2rootR))
        dt = (2r * K + 2r * (r^2 + a^2)) / D - (r^2 + a^2) * K * Dprime / D^2
        dphi = 2r * a / D - a * K * Dprime / D^2
        return eps(T) * abs(r) .* (abs(dt), abs(dphi), 2abs(r))
    end
    Δ = kerr_delta(a, r); Δprime = 2 * (r - 1)
    dt = (2r * P + (r^2 + a^2) * 2E * r) / Δ - (r^2 + a^2) * P * Δprime / Δ^2
    dphi = 2E * a / Δ - a * _rc_azimuth_numerator(c,r) * Δprime / Δ^2
    return eps(T) * abs(r) .* (abs(dt), abs(dphi), 2abs(r))
end

# ---- the end at infinity (header, :infinity) ------------------------------------------------
struct _InfinityTail{T}
    rs::T                                        # split radius (q = 1)
    W::NTuple{5,T}                               # c₄, β, γ, δ, ε of W(u)
    F3::NTuple{3,T}; F1::NTuple{3,T}             # coefficients of q⁻³, q⁻¹ (t, φ, τ)
    σ::T
end
_no_tail(::Type{T}) where {T} = (o = zero(T);
    _InfinityTail{T}(T(NaN), (o, o, o, o, o), (o, o, o), (o, o, o), one(T)))

function _InfinityTail(c::_RadialConstants{T}, rs, σ) where {T}
    c0, c1, c2, c3, c4 = c.R
    E = c.E
    o = zero(T)
    return _InfinityTail{T}(rs, (c4, c3 / rs, c2 / rs^2, c1 / rs^3, c0 / rs^4),
        (-2E * rs / σ, o, -2rs / σ), (-4E / σ, o, o), σ)
end

# artanh(x) given x and 1 − x (1 − x carries the digits when x → 1)
_atanh_stable(x, omx) = x < 0.5 ? atanh(x) : 0.5 * log((2 - omx) / omx)
# (artanh x − x)/x³ = 1/3 + x²/5 + x⁴/7 + … (x < 0.1: each term below 1/100 of the previous;
# 12 terms in Float64, proportionally more digits in T)
function _atanh_g1(x, omx)
    x < 0.1 || return (_atanh_stable(x, omx) - x) / x^3
    s = zero(x); p = one(x)
    for k in 1:_nterms(typeof(x), 12)
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
    I2 = -(iszero(κ) ? one(x) : _atanh_stable(x, omx) / x) / S # ∫ dq/(q S)
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
        iszero(c4) || return (zero(c4), zero(c4), zero(c4))
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
mutable struct _RadialSegment{T}
    kind::Symbol
    lo::T; hi::T                             # λ range
    σ::T                                     # sign of dr/dλ
    anchor_λ::T                              # λ where the segment's value equals `base`
    prim::ChebPieces{T}
    prim_anchor::NTuple{3,T}
    base::NTuple{3,T}
    rd::T                                    # :asymptote only
    multiplicity::Int                        # :asymptote only: 2, or 3 (variable s)
    rates_d::NTuple{3,T}
    tail::_InfinityTail{T}                   # :infinity only
    rstar_anchor::T; azimuth_anchor::T       # :horizon only
    hs::T                                    # :horizon only: sign P(r_h)
end

_needs_r(s::_RadialSegment) = s.kind !== :plain

function _segment_value(s::_RadialSegment{T}, c::_RadialConstants, λ, r, r_of=nothing) where {T}
    kind = s.kind
    if kind === :plain
        return s.base .+ (_eval3(s.prim, T(λ)) .- s.prim_anchor)
    elseif kind === :horizon
        F = _eval3(s.prim, T(λ)) .- s.prim_anchor
        w = s.σ * s.hs
        return (s.base[1] + F[1] + w * (_engine_rstar(c,r_of,λ,r) - s.rstar_anchor),
                s.base[2] + F[2] + w * (_engine_azimuth(c,r_of,λ,r) - s.azimuth_anchor),
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

function _tail_value(s::_RadialSegment{T}, r) where {T}
    q = min(T(_tail_q(s.tail, r)), one(T))
    return _tail_principal(s.tail, q) .+ _eval3(s.prim, q)
end

# the horizon-regular parts: t − sσ r*, φ − sσ φ_H (finite at the horizon this segment meets)
function _segment_regular(s::_RadialSegment{T}, λ) where {T}
    F = _eval3(s.prim, T(λ)) .- s.prim_anchor
    w = s.σ * s.hs
    return (s.base[1] + F[1] - w * s.rstar_anchor, s.base[2] + F[2] - w * s.azimuth_anchor)
end

function _build_segment(kind, c::_RadialConstants{T}, r_of, lo, hi, σ; rd=T(NaN),
        multiplicity=2, rh=c.rplus) where {T}
    a = c.a; E = c.E
    o = zero(T)
    tol = _radial_engine_tol(T)
    floors = (o, abs(a) * (1 + abs(E)), o)
    rates_d = (o, o, o); tail = _no_tail(T); hs = one(T)
    if kind === :plain
        f = λ -> _engine_rates(:plain,c,r_of,λ,hs)
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=tol, absfloor=floors,
            abserr=λ -> _engine_rounding(:plain,c,r_of,λ,hs)))
    elseif kind === :horizon
        hs = _rc_momentum(c,rh) >= 0 ? one(T) : -one(T)
        state = _radial_state(r_of,(lo+hi)/2)
        state === nothing || (hs = T(sign(state.chart.momentum)))
        f = λ -> _engine_rates(:horizon,c,r_of,λ,hs)
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=tol, absfloor=floors,
            abserr=λ -> _engine_rounding(:horizon,c,r_of,λ,hs)))
    elseif kind === :infinity
        # variable q = √(r_s/r) from the finite end (q = 1) to infinity (q = 0)
        tail = _InfinityTail(c, r_of(σ > 0 ? lo : hi), σ)
        prim = chebintegrate(chebfit(q -> _tail_rest(c, tail, q), o, one(T); ncomp=3,
            tol=tol, absfloor=floors))
    elseif kind === :asymptote
        aa, EE, LL = c.a, c.E, c.L
        num_t = [EE * aa^4 - aa * LL * aa^2, o, 2EE * aa^2 - aa * LL, o, EE]
        num_φ = aa .* [2EE * o - aa * LL, 2EE]             # a(2E r − aL)
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
            prim = chebintegrate(chebfit(g, o, sqrt(abs(rfar - rd)); ncomp=3,
                tol=tol, absfloor=floors))
        else
            f = function (r)
                w = σ * sgn / sqrt(max(_horner(R2, r), floatmin(T)))
                Δr = _horner(Δc, r)
                return (w * _horner(Qt, r) / (Δr * Δd), w * _horner(Qφ, r) / (Δr * Δd),
                    w * (r + rd))
            end
            prim = chebintegrate(chebfit(f, min(rd, rfar), max(rd, rfar); ncomp=3,
                tol=tol, absfloor=floors))
        end
    else
        error("unknown radial segment kind $kind")
    end
    return _RadialSegment{T}(kind, lo, hi, σ, T(NaN), prim, (o, o, o),
        (o, o, o), rd, multiplicity, rates_d, tail, T(NaN), T(NaN), hs)
end

function _anchor!(s::_RadialSegment{T}, c, r_of, λ, base) where {T}
    r = _needs_r(s) ? r_of(λ) : T(NaN)
    s.anchor_λ = λ
    s.prim_anchor = s.kind === :infinity ? _tail_value(s, r) :
        _eval3(s.prim, s.kind === :asymptote ? _asymptote_x(s, r) : λ)
    s.base = base
    if s.kind === :horizon
        s.rstar_anchor = _engine_rstar(c,r_of,λ,r)
        s.azimuth_anchor = _engine_azimuth(c,r_of,λ,r)
    end
    return s
end

# λ on a monotone leg where r(λ) = target: bisection between an inner point and the far
# end (an open endpoint is not evaluated; an infinite one is bracketed by doubling).
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
struct RadialEngine{T,F,P}
    c::_RadialConstants{T,P}
    r_of::F
    segments::Vector{_RadialSegment{T}}
    periodic::Bool
    period::T
    totals::NTuple{3,T}
    shift::NTuple{3,T}
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
    T = _float_type(a, E, L, Q)
    t = _radial_engine_tables(T(a), T(E), T(L), T(Q), _ErasedFunction{T}(r_of),
        _ErasedFunction{T}(potential), (T(domain[1]), T(domain[2])), (ends[1], ends[2]),
        turn === nothing ? nothing : T(turn), T(σ), T(rd), Int(multiplicity),
        T(λ_ref), period === nothing ? nothing : T(period))
    c = _rc(a, E, L, Q, potential)
    return RadialEngine{T,typeof(r_of),typeof(potential)}(c, r_of, t.segments, t.periodic,
        t.period, t.totals, t.shift)
end

# A radius or potential function behind a field of abstract type. The engine's tables are built
# through it, so the table construction (segment search, Chebyshev fits, anchoring) is compiled
# once for every radial model (per floating-point type) instead of once per model, at the price
# of a dynamic call per sample; the engine evaluates the tables with the concrete functions.
struct _ErasedFunction{T}
    f::Any
end
(e::_ErasedFunction{T})(x) where {T} = convert(T, e.f(x))::T
_radial_state(e::_ErasedFunction,lambda) = _radial_state(e.f,lambda)

function _radial_engine_tables(a::T, E::T, L::T, Q::T, r_of::_ErasedFunction{T},
        potential::_ErasedFunction{T}, domain, ends, turn, σ, rd, multiplicity, λ_ref,
        period) where {T}
    c = _rc(a, E, L, Q, potential)
    zero3 = (zero(T), zero(T), zero(T))
    engine(segments, periodic, period, totals, shift) =
        RadialEngine{T,typeof(r_of),typeof(potential)}(c, r_of, segments, periodic, period,
            totals, shift)
    if period !== nothing && turn !== nothing
        e = _radial_engine_tables(a, E, L, Q, r_of, potential, domain, (:turning, :turning), turn,
            σ, T(NaN), 2, domain[1], nothing)
        lo = domain[1]
        totals = _radial_eval_raw(e, lo + period)
        e = engine(e.segments, true, period, totals, zero3)
        v = _radial_eval_raw(e, λ_ref)
        return engine(e.segments, true, period, totals, v)
    elseif period !== nothing
        lo = domain[1]
        s = _build_segment(:plain, c, r_of, lo, lo + period, one(T))
        _anchor!(s, c, r_of, lo, zero3)
        totals = ntuple(k -> s.prim(lo + period, k) - s.prim_anchor[k], 3)
        e = engine([s], true, period, totals, zero3)
        v = _radial_eval_raw(e, λ_ref)
        return engine([s], true, period, totals, v)
    end
    λa, λb = domain
    legs = turn === nothing ? [(λa, λb, ends[1], ends[2], σ)] :
        [(λa, turn, ends[1], :turning, -σ), (turn, λb, :turning, ends[2], σ)]
    segs = _RadialSegment{T}[]
    for (lo, hi, klo, khi, sg) in legs
        append!(segs, _leg_segments(c, r_of, lo, hi, klo, khi, sg, rd, multiplicity, λ_ref))
    end
    sort!(segs; by=s -> s.lo)
    # anchor: the segment holding λ_ref is zero there; neighbours continue from it
    i0 = findfirst(s -> s.lo <= λ_ref <= s.hi, segs)
    i0 === nothing && error("radial engine: reference λ = $λ_ref outside the domain.")
    _anchor!(segs[i0], c, r_of, λ_ref, zero3)
    for j in i0+1:length(segs)
        b = segs[j].lo
        _anchor!(segs[j], c, r_of, b, _segment_value(segs[j-1], c, b, r_of(b),r_of))
    end
    for j in i0-1:-1:1
        b = segs[j].hi
        _anchor!(segs[j], c, r_of, b, _segment_value(segs[j+1], c, b, r_of(b),r_of))
    end
    return engine(segs, false, zero(T), zero3, zero3)
end

_end_radius(c, r_of, λ, kind) = kind === :horizon ? c.rplus : r_of(λ)

function _leg_segments(c::_RadialConstants{T}, r_of, lo, hi, klo, khi, σ, rd, multiplicity,
        λ_ref) where {T}
    segs = _RadialSegment{T}[]
    open_end(k) = k === :infinity || k === :asymptote      # ends that get their own segment
    # finite inner reference radius of the leg
    rin_lo = open_end(klo) ? T(NaN) : _end_radius(c, r_of, lo, klo)
    rin_hi = open_end(khi) ? T(NaN) : _end_radius(c, r_of, hi, khi)
    # A closed endpoint is a valid finite bracket, even on a very short Mino interval.
    inner = !open_end(klo) ? lo : !open_end(khi) ? hi : λ_ref
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
    crossed = open_end(klo) || open_end(khi) || hor_lo || hor_hi ? T[] :
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
            state = _radial_state(r_of,lambda_t)
            logspan = state === nothing ? log(rt/rh) : log1p(state.gap/rh)
            count=max(2,ceil(Int,logspan/log(2.0)))
            cuts=T[lambda_h]
            for j in 1:count-1
                if state === nothing
                    target=rh*exp(logspan*j/count)
                    push!(cuts,_leg_lambda(r_of,lambda_h,lambda_t,target))
                else
                    target=rh*expm1(logspan*j/count)
                    push!(cuts,_leg_gap_lambda(r_of,lambda_h,lambda_t,target))
                end
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

function _leg_gap_lambda(r_of,left,right,target)
    lo,hi=minmax(left,right)
    initial=_radial_state(r_of,lo).gap-target
    for _ in 1:200
        mid=(lo+hi)/2
        (mid==lo || mid==hi) && break
        if (_radial_state(r_of,mid).gap-target)*initial>0
            lo=mid
        else
            hi=mid
        end
    end
    return (lo+hi)/2
end

# horizons strictly between two radii, in the order met going from `ra` to `rb` (r₋ only for
# a ≠ 0, where it is a pole of the rates)
function _crossed_horizons(c::_RadialConstants{T}, ra, rb) where {T}
    lo, hi = minmax(ra, rb)
    h = kerr_horizons(c.a)
    hs = T[x for x in (h.rplus, h.rminus) if lo < x < hi && !(x == h.rminus && iszero(c.a))]
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
        return _segment_value(s, e.c, μ, _needs_r(s) ? e.r_of(μ) : NaN,e.r_of) .+ n .* e.totals
    end
    s = _find_segment(e, λ)
    return _segment_value(s, e.c, λ, _needs_r(s) ? e.r_of(λ) : NaN,e.r_of)
end
function _radial_eval_raw(e::RadialEngine, λ, r)
    if e.periodic
        length(e.segments) == 1 && return _radial_eval_raw(e, λ)
        n = floor((λ - e.segments[1].lo) / e.period)
        μ = λ - n * e.period
        return _segment_value(_find_segment(e, μ), e.c, μ, r,e.r_of) .+ n .* e.totals
    end
    return _segment_value(_find_segment(e, λ), e.c, λ, r,e.r_of)
end

"""(t_r, φ_r, τ_r) at λ (zero at the engine's reference λ); pass r if already known."""
_radial_eval(e::RadialEngine, λ) = _radial_eval_raw(e, λ) .- e.shift
_radial_eval(e::RadialEngine, λ, r) = _radial_eval_raw(e, λ, r) .- e.shift

function _increment_radii(e,lo,hi,left,right,radii)
    rlo=radii!==nothing && lo==left ? radii[1] : e.r_of(lo)
    rhi=radii!==nothing && hi==right ? radii[2] : e.r_of(hi)
    return rlo,rhi
end
function _empty_increment_segment(s,lo,hi,left,right,radii)
    hi<lo && return true
    hi>lo && return false
    return !(s.kind===:infinity && radii!==nothing && left==right &&
        radii[1]!=radii[2] && s.lo<=left<=s.hi)
end

function _radial_regular_interval(e::RadialEngine{T}, left, right, σr,radii=nothing) where {T}
    left==right && (radii===nothing || radii[1]==radii[2]) && return (zero(T),zero(T))
    right<left && return .-_radial_regular_interval(e,right,left,σr,
        radii===nothing ? nothing : reverse(radii))
    value = (zero(T), zero(T))
    for s in e.segments
        lo = max(left, s.lo); hi = min(right, s.hi)
        _empty_increment_segment(s,lo,hi,left,right,radii) && continue
        coefficient = -σr
        if s.kind === :plain || s.kind === :horizon
            part = ntuple(k -> _cheb_increment(s.prim, lo, hi, k), 2)
            s.kind === :horizon && (coefficient += s.σ * s.hs)
        elseif s.kind === :asymptote
            rlo,rhi=_increment_radii(e,lo,hi,left,right,radii)
            xlo = _asymptote_x(s, rlo); xhi = _asymptote_x(s, rhi)
            part = ntuple(k -> _cheb_increment(s.prim, xlo, xhi, k) +
                s.rates_d[k] * (hi - lo), 2)
        else
            rlo,rhi=_increment_radii(e,lo,hi,left,right,radii)
            qlo = min(T(_tail_q(s.tail, rlo)), one(T))
            qhi = min(T(_tail_q(s.tail, rhi)), one(T))
            p_lo = _tail_principal(s.tail, qlo); p_hi = _tail_principal(s.tail, qhi)
            part = ntuple(k -> p_hi[k] - p_lo[k] +
                _cheb_increment(s.prim, qlo, qhi, k), 2)
        end
        # The matching horizon kernel already is the regular rate, including at the horizon.
        if !iszero(coefficient)
            rlo,rhi=_increment_radii(e,lo,hi,left,right,radii)
            geometry=s.kind===:infinity ?
                (_rstar_all(e.c.a,rhi)-_rstar_all(e.c.a,rlo),
                    _horizon_azimuth(e.c.a,rhi)-_horizon_azimuth(e.c.a,rlo)) :
                (_engine_rstar(e.c,e.r_of,hi,rhi)-_engine_rstar(e.c,e.r_of,lo,rlo),
                    _engine_azimuth(e.c,e.r_of,hi,rhi)-_engine_azimuth(e.c,e.r_of,lo,rlo))
            part = part .+ coefficient .* geometry
        end
        value = value .+ part
    end
    return value
end

function _radial_regular_increment(e::RadialEngine, left, right, σr,radii=nothing)
    right<left && return .-_radial_regular_increment(e,right,left,σr,
        radii===nothing ? nothing : reverse(radii))
    e.periodic || return _radial_regular_interval(e,left,right,σr,radii)
    lo = e.segments[1].lo
    nl = floor((left - lo) / e.period); nr = floor((right - lo) / e.period)
    nl == nr && return _radial_regular_interval(e,left-nl*e.period,right-nr*e.period,σr,radii)
    return (nr - nl - 1) .* (e.totals[1], e.totals[2]) .+
        _radial_regular_interval(e, left - nl * e.period, lo + e.period, σr) .+
        _radial_regular_interval(e, lo, right - nr * e.period, σr)
end

function _radial_proper_interval(e::RadialEngine{T},left,right,radii=nothing) where {T}
    left==right && (radii===nothing || radii[1]==radii[2]) && return zero(T)
    right<left && return -_radial_proper_interval(e,right,left,
        radii===nothing ? nothing : reverse(radii))
    value=zero(T)
    for s in e.segments
        lo=max(left,s.lo); hi=min(right,s.hi)
        _empty_increment_segment(s,lo,hi,left,right,radii) && continue
        if s.kind===:plain || s.kind===:horizon
            value+=_cheb_increment(s.prim,lo,hi,3)
        elseif s.kind===:asymptote
            rlo,rhi=_increment_radii(e,lo,hi,left,right,radii)
            xlo=_asymptote_x(s,rlo); xhi=_asymptote_x(s,rhi)
            value+=_cheb_increment(s.prim,xlo,xhi,3)+s.rates_d[3]*(hi-lo)
        else
            rlo,rhi=_increment_radii(e,lo,hi,left,right,radii)
            qlo=min(T(_tail_q(s.tail,rlo)),one(T))
            qhi=min(T(_tail_q(s.tail,rhi)),one(T))
            value+=_tail_principal(s.tail,qhi)[3]-_tail_principal(s.tail,qlo)[3]+
                _cheb_increment(s.prim,qlo,qhi,3)
        end
    end
    return value
end

function _radial_proper_increment(e::RadialEngine,left,right,radii=nothing)
    right<left && return -_radial_proper_increment(e,right,left,
        radii===nothing ? nothing : reverse(radii))
    e.periodic || return _radial_proper_interval(e,left,right,radii)
    lo=e.segments[1].lo
    nl=floor((left-lo)/e.period); nr=floor((right-lo)/e.period)
    nl==nr && return _radial_proper_interval(e,left-nl*e.period,right-nr*e.period,radii)
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
    n = e.periodic && length(e.segments) > 1 ? floor((λ - e.segments[1].lo) / e.period) :
        zero(e.period)
    μ = λ - n * e.period
    s = _find_segment(e, μ)
    if s.kind === :horizon && s.σ * s.hs == σ_wanted
        v = _segment_regular(s, μ)
        return (v[1] + n * e.totals[1] - e.shift[1], v[2] + n * e.totals[2] - e.shift[2])
    end
    r = e.r_of(λ)
    t, φ, _ = _radial_eval(e, λ, r)
    return (t-σ_wanted*_engine_rstar(e.c,e.r_of,λ,r),
        φ-σ_wanted*_engine_azimuth(e.c,e.r_of,λ,r))
end

_spectral_summary(e::RadialEngine) = _spectral_summary(s.prim for s in e.segments)

"""Mean radial rates over one period (periodic engines)."""
_radial_mean_rates(e::RadialEngine) = e.totals ./ e.period

"""
    EngineCoordinates

The coordinates of one member: radial engine part + polar primitive (`polar_primitive(λ)`
gives the polar (t, φ, τ) increments). One object holds, once, the radius function, the polar
primitive, the potential, the engine's build parameters, the origin of the regular chart and
the lazily built radial engine (built on the first t/φ/τ call and published atomically, so
concurrent first calls are safe), so that the closures of a member capture one small object
instead of trees of closures. Evaluate it with
`_coords_t`, `_coords_phi`, `_coords_tau`, `_coords_tphitau`, `_coords_v`, `_coords_psi`,
`_coords_radial`, `_coords_spectral`; `_regular_chart(c, σ, λ_ref)` gives further regular charts.
"""
# the lazily built radial engine of one member
mutable struct _EngineCache{E}
    @atomic engine::Union{Nothing,E}
end

struct EngineCoordinates{T,F,P,V}
    a::T; E::T; L::T; Q::T
    r_of::F
    polar_primitive::P
    potential::V
    domain::Tuple{T,T}
    ends::Tuple{Symbol,Symbol}
    turn::Union{Nothing,T}
    σ::T
    rd::T
    multiplicity::Int
    period::Union{Nothing,T}
    λ_bl::T                         # t = φ = 0 here
    pb::NTuple{3,T}                 # polar primitive at λ_bl
    λ_tau::T                        # τ = 0 here
    pτ::T                           # polar τ primitive at λ_tau
    λ_regular::Union{Nothing,T}     # v = ψ = 0 here (nothing: no regular chart)
    σ_regular::T
    cache::_EngineCache{RadialEngine{T,F,V}}
    regular_origin::NTuple{2,T}     # polar (t, φ) primitive at λ_regular
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
    T = _float_type(a, E, L, Q)
    return EngineCoordinates{T,F,P,V}(a, E, L, Q, r_of, polar_primitive, potential,
        (T(domain[1]), T(domain[2])), (ends[1], ends[2]),
        turn === nothing ? nothing : T(turn), σ, rd, multiplicity,
        period === nothing ? nothing : T(period), λ_bl, T.(polar_primitive(λ_bl)), λ_tau,
        polar_primitive(λ_tau)[3], λ_regular === nothing ? nothing : T(λ_regular),
        σ_regular, _EngineCache{RadialEngine{T,F,V}}(nothing),
        λ_regular === nothing ? (T(NaN), T(NaN)) :
            T.(_polar_origin(polar_primitive, λ_regular)))
end

function _coords_engine(c::EngineCoordinates)
    e = @atomic c.cache.engine
    e === nothing || return e
    e = _radial_engine(c.a, c.E, c.L, c.Q, c.r_of; potential=c.potential, domain=c.domain,
        ends=c.ends, turn=c.turn, σ=c.σ, rd=c.rd, multiplicity=c.multiplicity, λ_ref=c.λ_bl,
        period=c.period)
    @atomic c.cache.engine = e
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

# the polar (t, φ) primitive at the origin λr of a regular chart
_polar_origin(polar_primitive, λr) = (p0 = polar_primitive(λr); (p0[1], p0[2]))

# Definite regular-coordinate increment from λr, with the polar origin at λr.
function _coords_regular(c::EngineCoordinates, σr, λr, origin, λ)
    e = _coords_engine(c)
    v = _radial_regular_increment(e, λr, λ, σr)
    p = c.polar_primitive(λ)
    return (v[1] + p[1] - origin[1], v[2] + p[2] - origin[2])
end
function _coords_regular(c::EngineCoordinates, λ)
    c.λ_regular === nothing && error("No horizon-regular reference event for this member.")
    return _coords_regular(c, c.σ_regular, c.λ_regular, c.regular_origin, λ)
end
_coords_v(c::EngineCoordinates, λ) = _coords_regular(c, λ)[1]
_coords_psi(c::EngineCoordinates, λ) = _coords_regular(c, λ)[2]

"""A further regular chart (t − σr r*, φ − σr φ_H), zero at λr."""
function _regular_chart(c::EngineCoordinates, σr, λr)
    origin = _polar_origin(c.polar_primitive, λr)
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
    radial(k) = (r1, r2) -> σ * (k==3 ?
        _radial_proper_increment(_coords_engine(coords),λ_of(r1),λ_of(r2),(float(r1),float(r2))) :
        _radial_regular_increment(_coords_engine(coords),λ_of(r1),λ_of(r2),zero(coords.a),(float(r1),float(r2)))[k])
    chart(k) = (r1, r2) -> σ *
        _radial_regular_increment(_coords_engine(coords),λ_of(r1),λ_of(r2),σ_regular,(float(r1),float(r2)))[k]
    increments = (radial_mino_increment=(r1, r2) -> σ * (λ_of(r2) - λ_of(r1)),
        radial_time_increment=radial(1), radial_phi_increment=radial(2),
        radial_proper_increment=radial(3))
    return regular ? merge(increments,
        (radial_v_increment=chart(1), radial_psi_increment=chart(2))) : increments
end
