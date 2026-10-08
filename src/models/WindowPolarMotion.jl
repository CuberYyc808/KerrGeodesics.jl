# Polar motion of E ≥ 1 scatter and capture orbits: `_window_polar_motion` (the pendular
# polar-engine solution with the E = 1 phase convention) and the E > 1 polar roots.


# The polar phase at the reference radial event and the polar event at phase zero, taken
# from `_window_polar_motion` (northern turning point for E > 1, northward equator crossing
# at E = 1, `:not_applicable` for equatorial orbits).
_polar_phase_metadata(polar, polar_phase, reference_event) =
    (polar_phase=polar.inclined ? polar_phase : nothing, reference_radial_event=reference_event,
     phase_convention=polar.metadata.phase_convention)

# E > 1 polar roots: z_+ in (0,1] (turning point) and z_- < 0, with m = z_+/(z_+ - z_-),
# 1 - m = -z_-/(z_+ - z_-) and omega = sqrt(beta (z_+ - z_-)), beta = a^2(E^2 - 1) = -c;
# formed from c*u so that nothing is divided by the small beta.
function _hyperbolic_polar_roots(a, energy, lz, q)
    roots = _polar_quadratic_roots(a, energy, lz, q)
    roots.disc >= 0 || error("Generic inclined hyperbolic polar roots are not real.")
    pairs = ((roots.u_small, roots.cu_small), (roots.u_big, roots.cu_big))
    (zplus, cplus), (zminus, cminus) = pairs[1][1] >= pairs[2][1] ? pairs : reverse(pairs)
    # (z_+ = 1 on the axis, Lz = 0; Θ(1) = −Lz² puts z_+ ≤ 1, and when Lz² ≪ Q the rounded
    # root can land an ulp above 1; the polar engine forms 1 − z_+ separately with
    # _polar_one_minus_root)
    # sqrt(u) remains representable when the small root u=Q/(c*u_big) does not.
    beta = -roots.c
    positive_large = roots.cu_big > 0
    amplitude = positive_large ? sqrt(q)/sqrt(roots.cu_big) :
        sqrt(-roots.cu_big)/sqrt(beta)
    negative_amplitude = positive_large ? sqrt(roots.cu_big)/sqrt(beta) :
        sqrt(q)/sqrt(-roots.cu_big)
    0 < amplitude <= sqrt(1 + _polar_root_slack(_float_type(amplitude))) || error("Generic inclined hyperbolic polar root z_+ is outside (0,1].")
    zplus = min(zplus, 1.0)
    negative_amplitude > 0 || error("Generic inclined hyperbolic polar root z_- must be negative.")
    spread = cminus - cplus                  # beta (z_+ - z_-) = -c z_+ + c z_-
    return (zplus=zplus, zminus=zminus, m=-cplus / spread, m1=cminus / spread,
        omega=sqrt(spread), amplitude=min(amplitude,1.0),
        complementary_modulus=negative_amplitude*sqrt(beta)/sqrt(spread))
end

"""
    _window_polar_motion(a, energy, lz, q, polar_phase)

Polar part of an E ≥ 1 orbit as Mino-time closures. Equatorial orbits (`Q = 0`) stay at
θ = π/2. Inclined orbits use the pendular polar-engine solution with phase 0 at the northern
turning point; at E = 1 the phase is shifted so that z = √z₋ sin(phase + ωλ), phase 0 on the
equator moving north. `polar_phase` is the phase at λ = 0.
`rootdata` is the polar entry of the orbit's `Roots` record.
"""
function _window_polar_motion(a, energy, lz, q, polar_phase)
    if iszero(q)
        T = _float_type(a, energy, lz, q)
        o = zero(T)
        return (inclined=false, z=λ -> o, uz=λ -> o, theta=λ -> T(π) / 2,
                phi=λ -> lz * λ, t=λ -> (a * lz - a^2 * energy) * λ, tau=λ -> o,
                phidot=λ -> lz, tdot=λ -> a * lz - a^2 * energy, rootdata=nothing,
                primitive=λ -> ((a * lz - a^2 * energy) * λ, lz * λ, o),
                position=λ -> (o, o, one(T)),
                metadata=(sector=:equatorial, phase=o, phase_convention=:not_applicable))
    end
    # polar engine (PolarEngine.jl): z = √z₊ cn(phase + ωλ) for E > 1, phase 0 at the
    # northern turning point; at E = 1 (exactly, `kerr_energy_regime`) the phase is shifted
    # by K so that z = √z₋ sin(phase + ωλ), phase 0 on the equator moving north
    parabolic = kerr_energy_regime(energy) === :parabolic
    sol = _polar_solution(a, energy, lz, q, :pendular, float(polar_phase))
    if parabolic
        sol = _polar_solution(a, energy, lz, q, :pendular,
            float(polar_phase) - sol.metadata.period_u / 2)
    end
    md = sol.metadata
    # z² roots: zplus the turning point, zminus the negative root (−∞ when β = a²(E² − 1)
    # vanishes and Θ is linear in z²)
    beta_zero = iszero(a^2 * _e2m1(energy))
    roots = _polar_quadratic_roots(a, energy, lz, q)
    zminus, zplus = beta_zero ? (-Inf, roots.u_small) : minmax(roots.u_small, roots.u_big)
    rootdata = (zplus=zplus, zminus=zminus, modulus=md.modulus, omega=md.omega,
        beta_zero=beta_zero)
    position = sol.position
    primitive = sol.primitive
    zf(λ) = position(λ)[1]
    return (inclined=true, z=zf, uz=λ -> position(λ)[2],
            theta=λ -> acos(clamp(zf(λ), -1.0, 1.0)), phi=λ -> primitive(λ)[2],
            t=λ -> primitive(λ)[1], tau=λ -> primitive(λ)[3],
            phidot=λ -> lz / position(λ)[3],
            tdot=λ -> a * lz - a^2 * energy * position(λ)[3],
            rootdata=rootdata, primitive=primitive, position=position,
            metadata=merge(md, (phase=float(polar_phase),
                phase_convention=parabolic ?
                    :equator_northward_at_zero_phase : :northern_turning_point_at_zero_phase)))
end
