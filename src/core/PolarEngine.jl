# Polar motion. Every oscillating polar sector has z(λ) = s·√A·J(u|m), u = u0 + ωλ, with
# J one of cd, cn, dn (cd = cn/dn; cd with m = 0 is cos). u = 0 is always the turning point
# nearest the axis, so there — where the φ rate is steepest — u carries relative, not
# absolute, rounding. The polar rates depend on z² only,
#     dt/dλ ⊃ aLz − a²E(1 − z²),   dφ/dλ ⊃ Lz/(1 − z²),   dτ/dλ ⊃ a²z²,
# and z² is periodic in u with period 2K(m). The primitives are therefore stored once per
# period as Chebyshev series in u (Spectral.jl) instead of being assembled from incomplete
# elliptic integrals at every call. sin²θ = 1 − z² is formed as ε + A·(complement of J²),
# ε = 1 − A from `_polar_one_minus_root`, and handed to the velocities as well.
#
# Near the axis the φ rate Lz/(1 − z²) is a spike of height Lz/ε and width √ε: its area
# (≈ π per pass) cannot be resolved by sampling once ε approaches the rounding of u. With a
# variable y(u) that vanishes at the spike, y' = 1 there and 1 − z² = ε + B y²,
#     Lz/(1 − z²) = Lz y'/(ε + B y²) + Lz (1 − y')/(ε + B y²),
# the first term integrates to Lz atan(√(B/ε) y)/√(εB) and the second is bounded (1 − y' =
# y²(p + q y²)/(1 + y')), so only it is fitted:
#     cd: y = sd u, y' = cn/dn², B = A k'²,  p = 1 − 2m, q = m k'²   (k'² = 1 − m)
#     cn: y = sn u, y' = cn dn,   B = A
#     dn: y = sn u, y' = cn dn,   B = A m,   p = 1 + m,  q = −m  (cn and dn)
# with the spike at u = 0 in every case.

"""
    _polar_one_minus_root(a, energy, lz, q, u)

`1 − u` for a root `u` of p(u) = c u² − (q + lz² + c) u + q (`c = a²(1 − E²)`), without
cancellation. The numbers y = 1 − u of both roots solve c y² − g y − lz² = 0 with
g = c − q − lz² (Vieta with p(1) = −lz²), so they come from a stable quadratic formula;
the one belonging to `u` is returned.
"""
function _polar_one_minus_root(a, energy, lz, q, u)
    c = -a^2 * _e2m1(energy)
    g = c - q - lz^2
    Y = g + copysign(sqrt(g^2 + 4c * lz^2), g)
    iszero(Y) && return 1 - u
    y_near = -2lz^2 / Y                      # the root nearer 1
    iszero(c) && return y_near
    y_far = Y / (2c)
    return abs(1 - y_near - u) <= abs(1 - y_far - u) ? y_near : y_far
end

# the spike variable y, its derivative y' and B, p, q of the φ-rate split (header)
@inline function _polar_spike(kind::Symbol, jac, A, m)
    sn, cn, dn = jac
    if kind === :cd
        kp2 = 1 - m
        return (sn / dn, cn / dn^2), (A * kp2, 1 - 2m, m * kp2)
    end
    return (sn, cn * dn), (kind === :cn ? A : A * m, 1 + m, -m)
end

# ∫ dy/(ε + B y²), including B = 0
_polar_spike_primitive(y, ε, B) = B > 0 ? atan(sqrt(B / ε) * y) / sqrt(ε * B) : y / ε

# z²/A and (1 − z²) for J(u|m); `jac` = (sn, cn, dn)
@inline function _polar_z2(kind::Symbol, jac, A, one_minus_A, m)
    sn, cn, dn = jac
    if kind === :cd                        # 1 − cd² = k'² sd²
        return A * (cn / dn)^2, one_minus_A + A * (1 - m) * (sn / dn)^2
    elseif kind === :cn
        return A * cn^2, one_minus_A + A * sn^2
    else                                   # :dn,  dn² = 1 − m sn²
        return A * dn^2, one_minus_A + A * m * sn^2
    end
end

"""
(sn, cn, dn)(u | m) with u first reduced to [−K, K] (sn, cn have period 4K and change sign
under u → u ± 2K), so large arguments and m → 1 keep full accuracy.
"""
function _ellipj_reduced(u, m, K)
    y = u - 4K * round(u / (4K))                 # [−2K, 2K]
    flip = false
    if y > K
        y -= 2K; flip = true
    elseif y < -K
        y += 2K; flip = true
    end
    sn, cn, dn = Elliptic.ellipj(y, m)
    return flip ? (-sn, -cn, dn) : (sn, cn, dn)
end

@inline function _polar_j(kind::Symbol, jac, m)
    sn, cn, dn = jac
    kind === :cd && return cn / dn, -(1 - m) * sn / dn^2
    kind === :cn && return cn, -sn * dn
    return dn, -m * sn * cn
end

"""
    _equatorial_polar_solution(a, energy, lz)

z ≡ 0: constant polar rates (t: aLz − a²E, φ: Lz).
"""
function _equatorial_polar_solution(a, energy, lz)
    rt = a * lz - a^2 * energy
    rates_primitive(lambda) = (rt * lambda, lz * lambda, 0.0)
    position(lambda) = (0.0, 0.0, 1.0)
    formula(lambda) = (z=0.0, uz=0.0, sin2=1.0, theta=pi / 2, phi=lz * lambda, t=rt * lambda,
        tau=0.0)
    return (formula=formula, primitive=rates_primitive, position=position,
        metadata=(sector=:equatorial, phase=0.0, phase_convention=:not_applicable,
            mean_rates=(t=rt, phi=float(lz), tau=0.0)))
end

"""
    _elliptic_polar_solution(a, energy, lz, q; kind, A, one_minus_A, m, omega, u0,
                             sign=1.0, metadata)

Polar motion z = sign·√A·J(u0 + ωλ | m) with t, φ, τ primitives from one Chebyshev period.
Returns `(formula, primitive, position, metadata)`: `formula(λ)` gives
`(z, uz = dz/dλ, sin2 = 1 − z², theta, phi, t, tau)`, `primitive(λ)` the polar `(t, φ, τ)`
and `position(λ)` `(z, dz/dλ, 1 − z²)`; the primitives vanish at λ = 0.
"""
function _elliptic_polar_solution(a, energy, lz, q; kind::Symbol, A, one_minus_A, m,
        omega, u0, sign=1.0, metadata::NamedTuple)
    K = Elliptic.K(m)
    period = 2K
    ε = one_minus_A
    _, (B, ps, qs) = _polar_spike(kind, (0.0, 0.0, 0.0), A, m)
    # the closed-form part of the φ primitive, zero at u = 0
    spike_primitive(u) = iszero(lz) ? 0.0 :
        lz / omega * _polar_spike_primitive(_polar_spike(kind, _ellipj_reduced(u, m, K), A, m)[1][1], ε, B)
    # its values at the ends of [0, K] are exact (y = 0 and 1/k' for cd, 1 for cn, dn): an
    # evaluated y at u = K would carry the rounding of K
    y_ends = kind === :cd ? (0.0, 1 / sqrt(1 - m)) : (0.0, 1.0)
    S0 = iszero(lz) ? 0.0 : lz / omega * _polar_spike_primitive(y_ends[1], ε, B)
    spike_half = iszero(lz) ? 0.0 : lz / omega * _polar_spike_primitive(y_ends[2], ε, B) - S0
    spike(u) = spike_primitive(u) - S0
    # z² is even about u = 0 and about u = K: fit the rates on [0, K] only (where ellipj is
    # accurate even for m → 1) and unfold by symmetry. The φ component is the bounded rest
    # of Lz/(1 − z²) after the closed-form spike.
    rates = function (u)
        jac = _ellipj_reduced(u, m, K)
        z2, omz2 = _polar_z2(kind, jac, A, one_minus_A, m)
        (y, yp), _ = _polar_spike(kind, jac, A, m)
        rest = iszero(lz) ? 0.0 : lz * y^2 * (ps + qs * y^2) / ((1 + yp) * omz2 * omega)
        return ((a * lz - a^2 * energy * omz2) / omega, rest, a^2 * z2 / omega)
    end
    scale_t = abs(a * lz) + a^2 * abs(energy) + 1.0e-300
    # the φ rest only needs the accuracy of the whole φ rate (mean |spike| rate over [0, K])
    fit = chebfit(rates, 0.0, K; ncomp=3,
        absfloor=(scale_t / omega, abs(spike_half) / K + 1.0e-300, a^2 / omega + 1.0e-300))
    prim = chebintegrate(fit)
    half = (chebtotal(prim, 1), chebtotal(prim, 2) + spike_half, chebtotal(prim, 3))
    totals = 2 .* half
    F(y) = (P = _eval3(prim, y); (P[1], P[2] + spike(y), P[3]))
    function primitive(u)
        n = floor(u / period)
        y = u - n * period
        Fy = y <= K ? F(y) : 2 .* half .- F(period - y)
        return n .* totals .+ Fy
    end
    p0 = primitive(float(u0))
    amp = sign * sqrt(A)
    # (t, φ, τ) polar primitives, and (z, dz/dλ, sin²θ = 1 − z² without cancellation)
    rates_primitive(lambda) = primitive(u0 + omega * float(lambda)) .- p0
    function position(lambda)
        jac = _ellipj_reduced(u0 + omega * float(lambda), m, K)
        J, dJ = _polar_j(kind, jac, m)
        return (amp * J, amp * omega * dJ, _polar_z2(kind, jac, A, one_minus_A, m)[2])
    end
    formula = function (lambda)
        z, uz, s2 = position(lambda)
        tφτ = rates_primitive(lambda)
        return (z=z, uz=uz, sin2=s2, theta=acos(clamp(z, -1.0, 1.0)),
            phi=tφτ[2], t=tφτ[1], tau=tφτ[3])
    end
    return (formula=formula, primitive=rates_primitive, position=position,
        metadata=merge(metadata,
        (modulus=m, omega=omega, formula_kind=Symbol(:spectral_jacobi_, kind),
         spectral=_spectral_summary((prim,)),
         period_u=period, mean_rates=(t=totals[1] / period * omega,
             phi=totals[2] / period * omega, tau=totals[3] / period * omega))))
end

"""
The polar solution `polar` delayed by δ: its value at λ is the original one at λ − δ, and the
primitives vanish at λ = 0 again. (A member whose λ origin moves by δ keeps its worldline.)
"""
function _polar_delayed(polar, δ)
    prim = _polar_primitive(polar)
    pos = _polar_position(polar)
    p0 = prim(-δ)
    primitive(λ) = prim(λ - δ) .- p0
    position(λ) = pos(λ - δ)
    function formula(λ)
        f = polar.formula(λ - δ)
        return merge(f, (t=f.t - p0[1], phi=f.phi - p0[2], tau=f.tau - p0[3]))
    end
    return (formula=formula, primitive=primitive, position=position,
        metadata=merge(polar.metadata, (phase_lambda=float(δ),)))
end

"""(achieved, pieces) of a polar solution's table (closed-form sectors have none)."""
_polar_spectral(polar) = hasproperty(polar, :metadata) && haskey(polar.metadata, :spectral) ?
    polar.metadata.spectral : (achieved=0.0, pieces=0)

"""(t, φ, τ) polar primitive of any polar solution (closed-form sectors via `formula`)."""
_polar_primitive(polar) = hasproperty(polar, :primitive) ? polar.primitive :
    (λ -> (f = polar.formula(λ); (f.t, f.phi, f.tau)))

"""(z, dz/dλ, sin²θ) of any polar solution (closed-form sectors give `sin2` in `formula`)."""
_polar_position(polar) = hasproperty(polar, :position) ? polar.position :
    (λ -> (f = polar.formula(λ); (f.z, f.uz, f.sin2)))
