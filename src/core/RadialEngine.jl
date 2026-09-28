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
struct _RadialConstants
    a::Float64; E::Float64; L::Float64; Q::Float64
    rplus::Float64
    R::Vector{Float64}          # radial potential coefficients
end

_rc(a, E, L, Q) = _RadialConstants(a, E, L, Q, kerr_horizons(a).rplus,
    collect(Float64, kerr_radial_coefficients(a, E, L, Q)))

@inline _rc_sqrtR(c::_RadialConstants, r) = sqrt(max(_horner(c.R, r), 0.0))

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
        br = s.prim.breaks
        x = clamp(r, min(br[1], br[end]), max(br[1], br[end]))
        return s.base .+ (_eval3(s.prim, x) .- s.prim_anchor) .+ s.rates_d .* (λ - s.anchor_λ)
    end
    return s.base .+ (_tail_value(s, r) .- s.prim_anchor)          # :infinity
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

function _build_segment(kind, c::_RadialConstants, r_of, lo, hi, σ; rd=NaN, rh=c.rplus)
    a = c.a; E = c.E
    floors = (0.0, abs(a) * (1 + abs(E)), 0.0)
    rates_d = (0.0, 0.0, 0.0); tail = _NO_TAIL; hs = 1.0
    if kind === :plain
        f = λ -> _plain_rates(c, r_of(λ))
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=_RADIAL_ENGINE_TOL, absfloor=floors))
    elseif kind === :horizon
        hs = kerr_radial_momentum(a, E, c.L, rh) >= 0 ? 1.0 : -1.0
        f = λ -> _horizon_rates(c, r_of(λ), hs)
        prim = chebintegrate(chebfit(f, lo, hi; ncomp=3, tol=_RADIAL_ENGINE_TOL, absfloor=floors))
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
        f = function (r)
            w = σ * sgn / sqrt(max(_horner(R2, r), floatmin()))
            Δr = _horner(Δc, r)
            return (w * _horner(Qt, r) / (Δr * Δd), w * _horner(Qφ, r) / (Δr * Δd), w * (r + rd))
        end
        prim = chebintegrate(chebfit(f, min(rd, rfar), max(rd, rfar); ncomp=3,
            tol=_RADIAL_ENGINE_TOL, absfloor=floors))
    else
        error("unknown radial segment kind $kind")
    end
    return _RadialSegment(kind, lo, hi, σ, NaN, prim, (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0), rd, rates_d, tail, NaN, NaN, hs)
end

function _anchor!(s::_RadialSegment, c, r_of, λ, base)
    r = _needs_r(s) ? r_of(λ) : NaN
    s.anchor_λ = λ
    s.prim_anchor = s.kind === :infinity ? _tail_value(s, r) :
        _eval3(s.prim, s.kind === :asymptote ? r : λ)
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
struct RadialEngine{F}
    c::_RadialConstants
    r_of::F
    segments::Vector{_RadialSegment}
    periodic::Bool
    period::Float64
    totals::NTuple{3,Float64}
    shift::NTuple{3,Float64}
end

"""
    _radial_engine(a, E, Lz, Q, r_of; domain, ends, turn=nothing, σ, rd=NaN, λ_ref=0.0)

`ends = (kind_lo, kind_hi)` with kinds `:turning` (any finite regular point), `:horizon`, `:infinity`,
`:asymptote`; `turn` is an interior turning point (or `nothing`) and `σ` the sign of
dr/dλ just above `turn` (or on the whole domain); `rd` the double root of an asymptote.
For a periodic libration pass `period` and `domain = (λ_peri, λ_peri + period)`; with a
`turn` as well (ends are then turning points) the period is cut into its two legs, which
may cross the horizons.
"""
function _radial_engine(a, E, L, Q, r_of; domain, ends=(:turning, :turning), turn=nothing,
        σ=1.0, rd=NaN, λ_ref=0.0, period=nothing)
    c = _rc(a, E, L, Q)
    if period !== nothing && turn !== nothing
        e = _radial_engine(a, E, L, Q, r_of; domain=domain, ends=(:turning, :turning),
            turn=turn, σ=σ, λ_ref=domain[1])
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
        append!(segs, _leg_segments(c, r_of, lo, hi, klo, khi, sg, rd, λ_ref))
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

function _leg_segments(c, r_of, lo, hi, klo, khi, σ, rd, λ_ref)
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
    open_end(klo) && push!(segs, _build_segment(klo, c, r_of, lo, s_lo, σ; rd=rd))
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
            m = (s_lo + s_hi) / 2
            hor_lo ? (push!(segs, _build_segment(:horizon, c, r_of, s_lo, m, σ));
                      push!(segs, _build_segment(:plain, c, r_of, m, s_hi, σ))) :
                     (push!(segs, _build_segment(:plain, c, r_of, s_lo, m, σ));
                      push!(segs, _build_segment(:horizon, c, r_of, m, s_hi, σ)))
        else
            push!(segs, _build_segment(:horizon, c, r_of, s_lo, s_hi, σ))
        end
    else
        push!(segs, _build_segment(:plain, c, r_of, s_lo, s_hi, σ))
    end
    open_end(khi) && push!(segs, _build_segment(khi, c, r_of, s_hi, hi, σ; rd=rd))
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
    _engine_coordinates(a, E, Lz, Q, r_of, polar_primitive; domain, ends, turn=nothing,
                        σ=1.0, rd=NaN, λ_bl=0.0, λ_tau=0.0, λ_regular=nothing,
                        σ_regular=-1.0, period=nothing)

Coordinates of one member as radial engine part + polar primitive (`polar_primitive(λ)`
gives the polar (t, φ, τ) increments). `t`, `phi` vanish at `λ_bl`, `tau` at `λ_tau` (the
λ = 0 event of every member), and
the horizon-regular pair `v = t − σ_regular r*`, `psi = φ − σ_regular φ_H` at `λ_regular`
(the horizon event, where they stay finite). The engine is built on the first call.
"""
function _engine_coordinates(a, E, L, Q, r_of::F, polar_primitive::P; domain, ends,
        turn=nothing, σ=1.0, rd=NaN, λ_bl=0.0, λ_tau=0.0, λ_regular=nothing,
        σ_regular=-1.0, period=nothing) where {F,P}
    cache = Ref{Union{Nothing,RadialEngine{F}}}(nothing)
    function engine()
        e = cache[]
        e === nothing || return e
        e = _radial_engine(a, E, L, Q, r_of; domain=domain, ends=ends, turn=turn, σ=σ,
            rd=rd, λ_ref=float(λ_bl), period=period)
        cache[] = e
        return e
    end
    pb = polar_primitive(λ_bl)
    tφτ(λ) = _radial_eval(engine(), λ) .+ (polar_primitive(λ) .- pb)
    τ0 = Ref(NaN)
    function tau(λ)
        isnan(τ0[]) && (τ0[] = tφτ(λ_tau)[3])
        return tφτ(λ)[3] - τ0[]
    end
    # horizon-regular pair (t − σr r*, φ − σr φ_H) + polar part, zero at λr
    function regular_chart(σr, λr)
        ref = Ref((NaN, NaN))
        return function (λ)
            e = engine()
            if isnan(ref[][1])
                v0 = _radial_regular(e, λr, σr)
                p0 = polar_primitive(λr)
                ref[] = (v0[1] + p0[1], v0[2] + p0[2])
            end
            v = _radial_regular(e, λ, σr)
            p = polar_primitive(λ)
            return (v[1] + p[1] - ref[][1], v[2] + p[2] - ref[][2])
        end
    end
    regular = λ_regular === nothing ?
        (λ -> error("No horizon-regular reference event for this member.")) :
        regular_chart(σ_regular, λ_regular)
    return (t=λ -> tφτ(λ)[1], phi=λ -> tφτ(λ)[2], tau=tau,
        v=λ -> regular(λ)[1], psi=λ -> regular(λ)[2], tphitau=tφτ,
        radial=λ -> _radial_eval(engine(), λ), regular_chart=regular_chart, engine=engine,
        spectral=() -> _spectral_summary(engine()))
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
    engine = coords.engine                  # capture only the engine accessor (small types)
    # the known radius is passed on: next to the horizon t, φ carry log|r − r₊|
    radial(k) = (r1, r2) -> σ * (_radial_eval(engine(), λ_of(r2), float(r2))[k] -
                                 _radial_eval(engine(), λ_of(r1), float(r1))[k])
    chart(k) = (r1, r2) -> σ * (_radial_regular(engine(), λ_of(r2), σ_regular)[k] -
                                _radial_regular(engine(), λ_of(r1), σ_regular)[k])
    increments = (radial_mino_increment=(r1, r2) -> σ * (λ_of(r2) - λ_of(r1)),
        radial_time_increment=radial(1), radial_phi_increment=radial(2),
        radial_proper_increment=radial(3))
    return regular ? merge(increments,
        (radial_v_increment=chart(1), radial_psi_increment=chart(2))) : increments
end
