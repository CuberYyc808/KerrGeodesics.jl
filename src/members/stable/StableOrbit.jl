# Class A (Stable) trajectories built directly from the constants of motion.
#
# r(λ) is the Jacobi-sn libration between the turning points r2 ≤ r ≤ r1 of the classified
# component, z(λ) the polar engine (PolarEngine.jl); t, φ and the Mino frequencies ϒt, ϒφ
# come from one Chebyshev period of the radial and polar rates. Nothing depends on inverting
# APEX (p, e, x) (e → 1, E → 1 stay well posed); with APEX input the component's roots are
# those of the turning-point geometry (`_apex_root_geometry`) when it is resolved.

"""Inner roots r3 ≥ r4 of R(r) once the turning points r1 ≥ r2 are known."""
function _class_a_inner_roots(a, energy, lz, q, r1, r2)
    c = [-a^2 * q, 2 * ((a * energy - lz)^2 + q), a^2 * (energy - 1) * (energy + 1) - lz^2 - q,
        2.0, (energy - 1) * (energy + 1)]
    c0, c1, c2 = _deflate_largest(_deflate_largest(c, r1), r2)
    disc = max(0.0, c1^2 - 4 * c2 * c0)
    big = -(c1 + copysign(sqrt(disc), c1)) / 2           # c2 r² + c1 r + c0, no cancellation
    rA = big / c2
    rB = iszero(big) ? zero(big) : c0 / big
    return max(rA, rB), min(rA, rB)
end

# the labels of kerr_geo_orbit_type_metadata: [family, shape, inclination]
_class_a_labels(case_id, q) = ["Stable", case_id === :A1 ? "Eccentric" : "Circular",
    iszero(q) ? "Equatorial" : "Inclined"]

"""
    _class_a_orbit(a, energy, lz, q, component; initPhases)

`KerrGeoStable` for a classified A1/A2 component. Phase conventions are those of
`kerr_geo_orbit`: at λ = 0 with zero phases the orbit is at periapsis and at its
northern polar turning point.
"""
function _class_a_orbit(a, energy, lz, q, component; initPhases=(0.0, 0.0, 0.0, 0.0))
    T = _float_type(a, energy, lz, q)
    a, energy, lz, q = T(a), T(energy), T(lz), T(q)
    case_id = component.CaseId
    r2 = float(component.LowerEndpoint.Radius)
    r1 = case_id === :A1 ? float(component.UpperEndpoint.Radius) : r2
    # APEX input: the inner roots of the turning-point geometry the component was classified with
    structure = component.Metadata.structure
    r3, r4 = haskey(structure, :apex_turning_points) ? structure.apex_turning_points.roots[3:4] :
        _class_a_inner_roots(a, energy, lz, q, r1, r2)
    potential = _radial_potential_from_roots(a, energy, lz, q, component.Metadata.structure)
    rc = _rc(a, energy, lz, q, potential)
    qt0, qr0, qθ0, qϕ0 = T.(initPhases)

    # ---- radial: r = r2 + (r1−r2)(r2−r3) sn² / ((r1−r3) cn² + (r2−r3) sn²), u = ωr λ
    # (A2: r1 = r2, k = 0 and ϒr is the epicyclic frequency)
    libration = case_id === :A1 ? _libration_model(energy,r1,r2,r3,r4) : nothing
    omega = case_id === :A1 ? libration.omega :
        sqrt(max(0.0,(1-energy)*(1+energy)*(r1-r3)*(r2-r4)))/2
    Kr = case_id === :A1 ? libration.K : _ellip_k(one(T))
    ϒr = π * omega / Kr
    radial_state(λ) = case_id === :A1 ?
        libration.state(λ) : (r2,zero(T))
    radial_r(λ) = radial_state(λ)[1]
    radial = case_id === :A1 ?
        _radial_engine(a, energy, lz, q, radial_r; potential=potential, domain=(zero(T), 2Kr / omega),
            period=2Kr / omega) :
        _plain_rates(rc, r2)

    # ---- polar: northern turning point at λ = 0
    polar = iszero(q) ? _equatorial_polar_solution(a, energy, lz) :
        _polar_solution(a, energy, lz, q, :pendular, zero(T))
    ϒθ = polar.metadata.sector === :equatorial ?
        sqrt(lz^2 + a^2 * (1 - energy) * (1 + energy)) :
        π * polar.metadata.omega / polar.metadata.period_u
    geometry = (a=a, energy=energy, lz=lz, q=q, case_id=case_id, roots=(r1, r2, r3, r4),
        rc=rc, ϒr=ϒr, ϒθ=ϒθ, phases=(qt0, qr0, qθ0, qϕ0))
    orbit, info = _class_a_assemble(geometry, radial_state, radial, polar.primitive,
        polar.position, values(polar.metadata.mean_rates))
    radial_tables = radial isa RadialEngine ? _spectral_summary(radial) : (achieved=0.0, pieces=0)
    return orbit, merge(info, (spectral=SpectralStatus(() -> (radial_tables, _polar_spectral(polar))),
        polar=polar.metadata))
end

# radial (t, φ, τ) primitive from periapsis: spectral for librations, linear for r = const
_class_a_radial_primitive(e::RadialEngine, λ) = _radial_eval(e, λ)
_class_a_radial_primitive(rates::NTuple{3,Real}, λ) = rates .* λ
_class_a_radial_mean(e::RadialEngine) = _radial_mean_rates(e)
_class_a_radial_mean(rates::NTuple{3,Real}) = rates

# (function barrier: every closure below captures concretely typed values)
function _class_a_assemble(g, radial_state::RS, radial::RE, polar_primitive::PP,
        polar_position::PZ, mean_θ) where {RS,RE,PP,PZ}
    (; a, energy, lz, q, case_id, rc, ϒr, ϒθ) = g
    r1, r2, r3, r4 = g.roots
    qt0, qr0, qθ0, qϕ0 = g.phases
    radial_primitive(λ) = _class_a_radial_primitive(radial, λ)
    radial_r(λ) = radial_state(λ)[1]
    mean_r = _class_a_radial_mean(radial)
    ϒt = mean_r[1] + mean_θ[1]
    ϒϕ = mean_r[2] + mean_θ[2]
    polar = (primitive=polar_primitive, position=polar_position)

    # ---- phases: λ offsets on the radial and polar clocks
    o = zero(ϒr)
    δr = ϒr > 0 ? qr0 / ϒr : o
    δθ = qθ0 / ϒθ
    R0 = radial_primitive(δr)
    Z0 = polar.primitive(δθ)
    t(λ) = qt0 + (radial_primitive(λ + δr)[1] - R0[1]) + (polar.primitive(λ + δθ)[1] - Z0[1])
    ϕ(λ) = qϕ0 + (radial_primitive(λ + δr)[2] - R0[2]) + (polar.primitive(λ + δθ)[2] - Z0[2])
    r(λ) = radial_r(λ + δr)
    θ(λ) = acos(clamp(polar.position(λ + δθ)[1], -1.0, 1.0))

    # ---- cross functions of the phases q_r = ϒr λ, q_θ = ϒθ λ (oscillating parts)
    Δtr(qr) = ϒr > 0 ? radial_primitive(qr / ϒr)[1] - mean_r[1] * qr / ϒr : o
    Δϕr(qr) = ϒr > 0 ? radial_primitive(qr / ϒr)[2] - mean_r[2] * qr / ϒr : o
    Δtθ(qθ) = polar.primitive(qθ / ϒθ)[1] - mean_θ[1] * qθ / ϒθ
    Δϕθ(qθ) = polar.primitive(qθ / ϒθ)[2] - mean_θ[2] * qθ / ϒθ
    dtr(qr) = ϒr > 0 ? (_plain_rates(rc, radial_r(qr / ϒr))[1] - mean_r[1]) / ϒr : o
    dϕr(qr) = ϒr > 0 ? (_plain_rates(rc, radial_r(qr / ϒr))[2] - mean_r[2]) / ϒr : o
    polar_rates(s2) = (a * lz - a^2 * energy * s2, iszero(lz) ? o : lz / s2)   # s2 = sin²θ
    dtθ(qθ) = (polar_rates(polar.position(qθ / ϒθ)[3])[1] - mean_θ[1]) / ϒθ
    dϕθ(qθ) = (polar_rates(polar.position(qθ / ϒθ)[3])[2] - mean_θ[2]) / ϒθ

    # ---- contravariant four-velocity u^μ = (dx^μ/dλ) / Σ
    function velocity(λ)
        rr, rdot = radial_state(λ + δr)
        z, uz, s2 = polar.position(λ + δθ)
        Σ = rr^2 + a^2 * z^2
        Tr, Φr, _ = _plain_rates(rc, rr)
        Tθ, Φθ = polar_rates(s2)
        return ((Tr + Tθ) / Σ, rdot / Σ, -uz / (Σ * sqrt(s2)), (Φr + Φθ) / Σ)
    end

    p = 2 * r1 * r2 / (r1 + r2)
    e = (r1 - r2) / (r1 + r2)
    x = _apex_x(a, energy, lz, q)
    orbit = KerrGeoStable(
        _class_a_labels(case_id, q),
        (a=a, p=p, e=e, x=x),
        (E=energy, Lz=lz, Q=q),
        "Mino",
        (t=t, r=r, θ=θ, ϕ=ϕ),
        (qt0=qt0, qr0=qr0, qθ0=qθ0, qϕ0=qϕ0),
        (ut=λ -> velocity(λ)[1], ur=λ -> velocity(λ)[2], uθ=λ -> velocity(λ)[3],
         uϕ=λ -> velocity(λ)[4]),
        (ϒt=ϒt, ϒr=ϒr, ϒθ=ϒθ, ϒϕ=ϒϕ),
        (Δtr=Δtr, Δtθ=Δtθ, Δϕr=Δϕr, Δϕθ=Δϕθ),
        (dtr=dtr, dtθ=dtθ, dϕr=dϕr, dϕθ=dϕθ),
    )
    # the member's own fields: z, dz/dλ and τ (zero at λ = 0) besides the orbit's t, r, θ, φ
    z(λ) = polar.position(λ + δθ)[1]
    τ(λ) = (radial_primitive(λ + δr)[3] - R0[3]) + (polar.primitive(λ + δθ)[3] - Z0[3])
    functions = (t=t, r=r, theta=θ, z=z, phi=ϕ, tau=τ, position=λ -> polar.position(λ + δθ),
        uz=λ -> polar.position(λ + δθ)[2],
        sin2=λ -> polar.position(λ + δθ)[3],
        sign_r=λ -> sign(radial_state(λ + δr)[2]))
    return orbit, (roots=g.roots, apex=(a=a, p=p, e=e, x=x), functions=functions)
end
