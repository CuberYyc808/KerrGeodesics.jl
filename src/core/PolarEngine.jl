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
#
# Lz = 0: the orbit passes over the axis, where φ is undefined. φ is the Lz → 0⁺ limit of the
# spike (x → 0⁺), a step of +π at every pass: Lz/(ω√(εB)) → 1 in every sector (Lz² = ε(Q − cA)/A
# at the turning root A), and atan(√(B/ε) y) → (π/2) sign(y).
#
# Close to the equator (small |Q| with Lz² < a²(E² − 1)) the parameter m approaches 1 and
# K ≈ ½ log(16/k'²) is fixed by k'² = 1 − m, which a rounded m no longer carries. Every
# sector therefore supplies m1 = k'² from its polar roots, and sn, cn, dn come from the
# Landen sequence started at √m1 (`_landen`, `_ellipj_reduced`).

"""
    _polar_one_minus_root(a, energy, lz, q, u)

`1 − u` for a root `u` of p(u) = c u² − (q + lz² + c) u + q (`c = a²(1 − E²)`), without
cancellation. The numbers y = 1 − u of both roots solve c y² − g y − lz² = 0 with
g = c − q − lz² (Vieta with p(1) = −lz²), so they come from a stable quadratic formula;
the one belonging to `u` is returned.
"""
function _polar_one_minus_root(a, energy, lz, q, u)
    roots = _polar_quadratic_roots(a, energy, lz, q)
    c = roots.c
    g = c - q - lz^2
    # The complementary quadratic has the same discriminant as the original one.
    square_root = abs(roots.cu_big-roots.cu_small)
    Y = g + copysign(square_root, g)
    iszero(Y) && return 1 - u
    y_near = -2lz^2 / Y                      # the root nearer 1
    iszero(c) && return y_near
    y_far = Y / (2c)
    return abs(1 - y_near - u) <= abs(1 - y_far - u) ? y_near : y_far
end

# the spike variable y, its derivative y' and B, p, q of the φ-rate split (header);
# `L` is the parameter record of `_landen`
@inline function _polar_spike(kind::Symbol, jac, A, L)
    sn, cn, dn = jac
    if kind === :cd
        return (sn / dn, cn / dn^2), (A * L.m1, 1 - 2L.m, L.m * L.m1)
    end
    return (sn, cn * dn), (kind === :cn ? A : A * L.m, 1 + L.m, -L.m)
end

# z²/A and (1 − z²) for J(u|m); `jac` = (sn, cn, dn)
@inline function _polar_z2(kind::Symbol, jac, A, one_minus_A, L)
    sn, cn, dn = jac
    if kind === :cd                        # 1 − cd² = k'² sd²
        return A * (cn / dn)^2, one_minus_A + A * L.m1 * (sn / dn)^2
    elseif kind === :cn
        return A * cn^2, one_minus_A + A * sn^2
    else                                   # :dn,  dn² = 1 − m sn²
        return A * dn^2, one_minus_A + A * L.m * sn^2
    end
end

@inline function _polar_j(kind::Symbol, jac, L)
    sn, cn, dn = jac
    kind === :cd && return cn / dn, -L.m1 * sn / dn^2
    kind === :cn && return cn, -sn * dn
    return dn, -L.m * sn * cn
end

function _polar_spike_increment(kind, L, left, delta, A, B, epsilon, lz_over_omega)
    right = left + delta
    jl = _ellipj_reduced(left, L); jr = _ellipj_reduced(right, L)
    if iszero(lz_over_omega)                # Lz → 0⁺: +π per pass over the axis
        yl = _polar_spike(kind, jl, A, L)[1][1]
        yr = _polar_spike(kind, jr, A, L)[1][1]
        return oftype(yl, π) / 2 * (sign(yr) - sign(yl))
    end
    mid = _ellipj_reduced(left + delta / 2, L)
    step = _ellipj_reduced(delta / 2, L)
    denominator = 1 - L.m * mid[1]^2 * step[1]^2
    difference = kind === :cd ?
        2step[1] * mid[2] * step[3] / (denominator * jl[3] * jr[3]) :
        2step[1] * mid[2] * mid[3] / denominator
    B > 0 || return lz_over_omega * difference / epsilon
    yl = _polar_spike(kind, jl, A, L)[1][1]
    yr = _polar_spike(kind, jr, A, L)[1][1]
    scale = sqrt(epsilon * B)
    angle = atan(scale * difference, epsilon + B * yl * yr)
    return lz_over_omega * angle / scale
end

function _polar_rates_increment(prim, totals, kind, L, initial, delta, A, B, epsilon,
        lz_over_omega)
    if iszero(delta)
        o = zero(_float_type(A, delta))
        return (o, o, o)
    end
    K = L.K; period = 2K
    cycles = round(delta / period)
    remainder = delta - cycles * period
    phase = initial - round(initial / period) * period
    value = cycles .* totals
    direction = remainder > 0 ? 1 : -1
    remaining = abs(remainder)
    # A residual interval spans at most two quarter-period boundaries.
    while remaining > 0
        quarter = direction > 0 ? floor(Int, phase / K) : ceil(Int, phase / K) - 1
        offset = phase - quarter * K
        orientation = iseven(quarter) ? 1 : -1
        local_phase = orientation > 0 ? offset : K - offset
        room = direction > 0 ? K - offset : offset
        width = min(remaining, max(room, 0.0))
        local_delta = orientation * direction * width
        part = ntuple(k -> _cheb_increment_delta(prim, local_phase, local_delta, k), 3)
        spike = _polar_spike_increment(kind, L, local_phase, local_delta, A, B, epsilon,
            lz_over_omega)
        value = value .+ orientation .* (part[1], part[2] + spike, part[3])
        remaining -= width
        remaining > 0 && (phase = direction > 0 ? (quarter + 1) * K : quarter * K)
    end
    return value
end

"""
    _equatorial_polar_solution(a, energy, lz)

z ≡ 0: constant polar rates (t: aLz − a²E, φ: Lz).
"""
function _equatorial_polar_solution(a, energy, lz)
    T = _float_type(a, energy, lz)
    rt = a * lz - a^2 * energy
    o = zero(T)
    rates_primitive(lambda) = (rt * lambda, lz * lambda, o)
    position(lambda) = (o, o, one(T))
    formula(lambda) = (z=o, uz=o, sin2=one(T), theta=T(π) / 2, phi=lz * lambda, t=rt * lambda,
        tau=o)
    return (formula=formula, primitive=rates_primitive, position=position,
        metadata=(sector=:equatorial, phase=o, phase_convention=:not_applicable,
            mean_rates=(t=rt, phi=T(lz), tau=o)))
end

"""
    _elliptic_polar_solution(a, energy, lz, q; kind, A, one_minus_A, m, m1, omega, u0,
                             sign=1.0, amplitude=sqrt(A),
                             complementary_modulus=sqrt(m1), metadata)

Polar motion z = sign·√A·J(u0 + ωλ | m) with t, φ, τ primitives from one Chebyshev period;
`m1` is 1 − m formed from the polar roots. `amplitude` carries sqrt(A) separately
when its square underflows. `complementary_modulus` similarly retains sqrt(m1).
Returns `(formula, primitive, position, metadata)`: `formula(λ)` gives
`(z, uz = dz/dλ, sin2 = 1 − z², theta, phi, t, tau)`, `primitive(λ)` the polar `(t, φ, τ)`
and `position(λ)` `(z, dz/dλ, 1 − z²)`; the primitives vanish at λ = 0.
"""
function _elliptic_polar_solution(a, energy, lz, q; kind::Symbol, A, one_minus_A, m, m1,
        omega, u0, sign=1.0, amplitude=sqrt(A), complementary_modulus=sqrt(m1),
        metadata::NamedTuple)
    T = _float_type(a, energy, lz, q, m, m1)
    L = _landen(T(m), T(m1); complementary_modulus=T(complementary_modulus))
    K = L.K
    period = 2K
    ε = one_minus_A
    _, (B, ps, qs) = _polar_spike(kind, (zero(T), zero(T), zero(T)), A, L)
    lz_over_omega = lz / omega
    spike_scale = !iszero(lz) && B > 0 ? sqrt(B / ε) : zero(T)
    spike_denom = !iszero(lz) && B > 0 ? sqrt(ε * B) : one(T)
    @inline spike_of_y(y) = iszero(lz) ? T(π) / 2 * Base.sign(y) : lz_over_omega *
        (B > 0 ? atan(spike_scale * y) / spike_denom : y / ε)
    # its values at the ends of [0, K] are exact (y = 0 and 1/k' for cd, 1 for cn, dn): an
    # evaluated y at u = K would carry the rounding of K
    y_ends = kind === :cd ? (zero(T), 1 / sqrt(T(m1))) : (zero(T), one(T))
    S0 = spike_of_y(y_ends[1])
    spike_half = spike_of_y(y_ends[2]) - S0
    # z² is even about u = 0 and about u = K: fit the rates on [0, K] only and unfold by
    # symmetry. The φ component is the bounded rest
    # of Lz/(1 − z²) after the closed-form spike.
    rates = function (u)
        jac = _ellipj_reduced(u, L)
        z2, omz2 = _polar_z2(kind, jac, A, one_minus_A, L)
        J = _polar_j(kind, jac, L)[1]
        (y, yp), _ = _polar_spike(kind, jac, A, L)
        rest = iszero(lz) ? zero(T) : lz * y^2 * (ps + qs * y^2) / ((1 + yp) * omz2 * omega)
        return ((a * lz - a^2 * energy * omz2) / omega, rest, J^2)
    end
    scale_t = abs(a * lz) + a^2 * abs(energy) + 1.0e-300
    # the φ rest only needs the accuracy of the whole φ rate (mean |spike| rate over [0, K])
    fit = chebfit(rates, zero(T), K; ncomp=3,
        absfloor=(scale_t / omega, abs(spike_half) / K + 1.0e-300, one(T)))
    prim = chebintegrate(fit)
    half = (chebtotal(prim, 1), chebtotal(prim, 2) + spike_half, chebtotal(prim, 3))
    totals = 2 .* half
    amp = sign * amplitude
    # Fit the dimensionless proper-time integrand before restoring its small amplitude.
    tau_scale = a^2 * amplitude * (amplitude / omega)
    # (t, φ, τ) polar primitives, and (z, dz/dλ, sin²θ = 1 − z² without cancellation)
    @inline function rates_primitive(lambda)
        value = _polar_rates_increment(prim, totals, kind, L, float(u0),
            omega * float(lambda), A, B, ε, lz_over_omega)
        return (value[1], value[2], tau_scale * value[3])
    end
    function position(lambda)
        jac = _ellipj_reduced(u0 + omega * float(lambda), L)
        J, dJ = _polar_j(kind, jac, L)
        return (amp * J, amp * omega * dJ, _polar_z2(kind, jac, A, one_minus_A, L)[2])
    end
    formula = function (lambda)
        z, uz, s2 = position(lambda)
        tφτ = rates_primitive(lambda)
        return (z=z, uz=uz, sin2=s2, theta=atan(sqrt(s2), z),
            phi=tφτ[2], t=tφτ[1], tau=tφτ[3])
    end
    return (formula=formula, primitive=rates_primitive, position=position,
        metadata=merge(metadata,
        (modulus=m, omega=omega, formula_kind=Symbol(:spectral_jacobi_, kind),
         spectral=_spectral_summary((prim,)),
         period_u=period, mean_rates=(t=totals[1] / period * omega,
             phi=totals[2] / period * omega, tau=tau_scale * totals[3] / period * omega))))
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
