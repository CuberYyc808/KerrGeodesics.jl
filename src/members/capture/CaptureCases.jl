# Class C (Capture) members C1–C12: from infinity into the future horizon. Radial models of C2
# and C4, the closed-form C1/C3 member, the shared model-based member (C2, C4, C5, C6–C12) and
# the exported constructors of the polar sectors.

const _NONGENERIC_POLAR_SECTORS = (
    :vortical, :constant_latitude, :equator_attractive, :axis_crossing)

function _c2_model(roots)
    x1, x2, x3 = roots
    leg = _three_real_leg((x1=x1, x2=x2, x3=x3))
    aroot, m, m1, scale = leg.A, leg.m, leg.m1, leg.scale
    phi_of_r(r) = _three_real_amplitude(leg, r)
    w(phi) = sqrt(max(1 - m * sin(phi)^2, 0.0))
    j2(phi) = -cos(phi) / sin(phi) * w(phi) + m * _ellip_d(phi, m1)       # F − E = m D
    function j4(phi)
        return ((2 + 2 * m) * j2(phi) - m * _ellip_f(phi, m1) -
            cos(phi) / sin(phi) * w(phi) / sin(phi)^2) / 3
    end
    function basis(r)
        phi = phi_of_r(r)
        f = _ellip_f(phi, m1)
        return (
            I0=-scale * f,
            I1=-scale * (x1 * f + aroot * j2(phi)),
            I2=-scale * (x1^2 * f + 2 * x1 * aroot * j2(phi) +
                aroot^2 * j4(phi)),
        )
    end
    pole(h, r) = -_three_real_pole(leg, h, phi_of_r(r))
    return (
        kind=:c2_parabolic_three_real_direct_capture,
        roots=(x1=x1, x2=x2, x3=x3),
        basis=basis,
        pole=pole,
        radius=δ -> _three_real_radius_from_infinity(leg, δ),
        mino=r -> _three_real_mino_from_infinity(leg, r),
        inward=true,
    )
end

# C4 has the Legendre form of D2 (members/scatter/FiniteWindow.jl) with (r_A, …, r_D) =
# (x1, …, x4): φ from its tangent and r from the Mino time δ after the infinity endpoint
# through `_d2_radius`. I0 is measured from infinity (Carlson's form, `_mino_to_infinity`),
# so the horizon-to-infinity time is not a difference of two F values.
function _c4_model(energy, roots)
    x1, x2, x3, x4 = roots
    leg = _four_real_leg(energy, (rA=x1, rB=x2, rC=x3, rD=x4))
    m1, scale = leg.m1, leg.prefactor
    # the characteristic of the leg's own factor 1/(r − x3) and its complement, n > 1
    n = (x4 - x1) / (x3 - x1)
    n1 = (x3 - x4) / (x3 - x1)
    lead = _e2m1(energy)
    function basis(r)
        phi = _four_real_amplitude(leg, r)
        f = _ellip_f(phi, m1)
        pin = _ellip_pi(phi, m1, n, n1)
        return (
            I0=-_mino_to_infinity(lead, roots, r),
            I1=scale * (x3 * f + (x4 - x3) * pin),
            I2=scale * (x3^2 * f + 2 * x3 * (x4 - x3) * pin +
                (x4 - x3)^2 * _ellip_pi2(phi, m1, n, n1)),
        )
    end
    function pole(h, r)
        phi = _four_real_amplitude(leg, r)
        nh, nh1 = _four_real_characteristic(leg, h)
        return scale * (_ellip_f(phi, m1) / (x3 - h) +
            (x3 - x4) * _ellip_pi(phi, m1, nh, nh1) / ((x4 - h) * (x3 - h)))
    end
    return (
        kind=:c4_hyperbolic_four_real_direct_capture,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        basis=basis,
        pole=pole,
        radius=δ -> _four_real_radius_from_infinity(leg, δ),
        mino=r -> _mino_to_infinity(lead, roots, r),
        inward=true,
    )
end

function _capture_analytic_model(case_id, energy, radii)
    case_id === :C2 && return _c2_model((radii[1], radii[2], radii[3]))
    case_id === :C4 && return _c4_model(
        energy, (radii[1], radii[2], radii[3], radii[4]))
    error("The analytic radial model covers C2 and C4, not $(case_id).")
end

# C1 (E = 1, one real root) and C3 (E > 1, two real roots): the closed-form radial models of
# the finite-window capture API. The pendular and equatorial sectors take the polar motion of
# `kerr_geo_capture`, the other sectors that of the polar engine.
function _closed_form_capture_component(a, energy, lz, q, component, polar_phase,
        polar_hemisphere, reference_radius)
    structure=component.Metadata.structure
    model = component.CaseId === :C1 ? _c1_radial_model(a, lz, q;structure=structure) :
        _c3_radial_model(a, energy, lz, q;structure=structure)
    model === nothing && error(component.CaseId === :C1 ?
        "The radial roots do not have the C1 structure (one real root below r₊ and a complex pair)." :
        "The radial roots do not have the C3 structure (two real roots below r₊ and a complex pair).")
    window = !(component.PolarSector in _NONGENERIC_POLAR_SECTORS)
    polar = window ? _window_polar_motion(a, energy, lz, q, polar_phase) :
        _polar_solution(a, energy, lz, q, component.PolarSector, polar_phase;
            hemisphere=polar_hemisphere)
    return _direct_component(a, energy, lz, q, component, model, polar, reference_radius)
end

# C1–C12: a radial model from infinity into the future horizon (λ = 0, where
# τ = v = ψ = 0); t and φ are zero at the reference radius (default: λ = λ_∞/2). All radial
# increments (t, φ, τ, v, ψ) come from the spectral engine.
function _direct_component(a, energy, lz, q, component, model, polar, reference_radius;
        roots=model.roots, status=(;))
    rplus = kerr_horizons(a).rplus
    λ_inf = -model.mino(rplus)
    λ_inf < 0 || error("The infinity endpoint must precede the future horizon.")
    check, check_bl = _infinity_to_horizon_checks(λ_inf)
    # r from the Mino time after the infinity endpoint, where the models keep their digits
    radius, mino = model.radius, model.mino
    r_of(λ) = λ == 0 ? rplus : radius(λ - λ_inf)
    function lambda_of_radius(r)
        r >= rplus || throw(DomainError(r, "The radius must lie outside the future horizon."))
        return λ_inf + mino(r)
    end
    λ_ref = reference_radius === nothing ? 0.5 * λ_inf : lambda_of_radius(float(reference_radius))
    λ_inf < λ_ref < 0 || error("The reference radius must be a finite exterior point.")
    potential = _radial_potential_from_roots(a, energy, lz, q, component.Metadata.structure)
    coords = _engine_coordinates(a, energy, lz, q, r_of, _polar_primitive(polar);
        potential=potential, domain=(λ_inf, 0.0), ends=(:infinity, :horizon), σ=-1.0, λ_bl=λ_ref,
        λ_regular=0.0, σ_regular=-1.0)
    engine = _radius_increments(coords, lambda_of_radius, -1.0)
    track = (r=λ -> r_of(check(λ)), check=check, check_bl=check_bl, coords=coords,
        sign_r=λ -> -1.0,
        domain=(mino=(λ_inf, 0.0), endpoint_closed=(false, true),
            endpoint_roles=(:past_infinity, :future_horizon), horizon_lambda=0.0),
        reference=(lambda0_event=:future_horizon, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=λ_ref, t_phi_zero_radius=r_of(λ_ref),
            tau_zero_event=:future_horizon, lambda_regular=0.0,
            polar_phase=polar.metadata.phase,
            polar_phase_convention=polar.metadata.phase_convention),
        trajectory=merge((lambda_of_radius=lambda_of_radius,), engine))
    return _engine_member(:capture, component.CaseId, a, energy, lz, q, polar, track,
        component, component.Metadata.structure, potential, kerr_geo_tier(component.CaseId),
        (model=roots,), merge((formula_family=component.FormulaFamily, formula_kind=model.kind),
            status), nothing)
end

# Mino-time checks of a member from infinity (λ_∞, open) into the horizon (λ = 0): (r, z, τ,
# v, ψ) on (λ_∞, 0], (t, φ) on (λ_∞, 0).
function _infinity_to_horizon_checks(λ_inf)
    check(λ) = (λ_inf < λ <= 0 || throw(DomainError(λ,
        "Mino time must lie in (λ_∞, 0] = ($(λ_inf), 0].")); float(λ))
    check_bl(λ) = (check(λ) < 0 || throw(DomainError(λ,
        "BL t and phi exclude the infinity and future-horizon endpoints.")); float(λ))
    return check, check_bl
end

"""
    kerr_geo_capture_component(a, E, Lz, Q; case_id=nothing, ...)

Classify (a, E, Lz, Q) and construct its Class C (Capture) component: from infinity into the
future horizon, with λ = 0 on the horizon (where τ = v = ψ = 0) and t, φ zero at
`reference_radius` (default `λ = λ_∞/2`). C1 and C3 use the closed-form radial formulas of
`kerr_geo_capture`; C2, C4 and C6–C12 their radial models (C6–C12: repeated roots strictly
below the outer horizon); C5 the four-complex-root model (`kerr_geo_capture_four_complex`).
Orbits from infinity that approach a Critical root asymptotically (K7, K10) are Critical
members, built by `kerr_geo_critical_component`.
"""
function kerr_geo_capture_component(a::Real, energy::Real, lz::Real, q::Real;
        case_id=nothing,
        polar_sector=nothing,
        polar_phase::Real=0.0,
        polar_hemisphere::Symbol=:north,
        reference_radius=nothing)
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=polar_sector)
    return _capture_component(a, energy, lz, q, classification, case_id;
        polar_sector=polar_sector, polar_phase=polar_phase, polar_hemisphere=polar_hemisphere,
        reference_radius=reference_radius)
end

# the Capture member `case_id` of classified constants
function _capture_component(a, energy, lz, q, classification, case_id; polar_sector=nothing,
        polar_phase=0.0, polar_hemisphere::Symbol=:north, reference_radius=nothing)
    component = _class_component(classification, :capture, case_id)
    phase = float(polar_phase)
    id = component.CaseId
    if id === :C5
        # C5 needs a polar sector even when the classifier leaves the choice to initial data
        classification = _with_polar_sector(classification,
            _resolved_polar_sector(classification, polar_sector, energy, lz, q))
        return _four_complex_component(a, energy, lz, q, classification,
            _class_component(classification, :capture, :C5); polar_hemisphere=polar_hemisphere,
            polar_phase=phase, reference_radius=reference_radius)
    end
    id in (:C1, :C3) && return _closed_form_capture_component(a, energy, lz, q, component,
        phase, polar_hemisphere, reference_radius)
    model = id in (:C2, :C4) ?
        _capture_analytic_model(id, energy, _root_radii(classification)) :
        interior_repeated_radial_model(id, a, energy, lz, q,
            classification.Status.root_structure)
    polar = _polar_solution(a, energy, lz, q, component.PolarSector, phase;
        hemisphere=polar_hemisphere)
    return _direct_component(a, energy, lz, q, component, model, polar, reference_radius)
end

function kerr_geo_capture_component(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_capture_component(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_capture_component(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_capture_component(a, constants...; kwargs...)
end

"""
    kerr_geo_capture_vortical(a, E, Lz, Q; kwargs...)

The Capture member of `(a, E, Lz, Q)` with vortical polar motion (E > 1, Q < 0: the orbit
stays in one hemisphere and never crosses the equator). Keywords as in
`kerr_geo_capture_component`.
"""
kerr_geo_capture_vortical(a::Real, energy::Real, lz::Real, q::Real; kwargs...) =
    kerr_geo_capture_component(a, energy, lz, q; polar_sector=:vortical, kwargs...)

"""
    kerr_geo_capture_constant_latitude(a, E, Lz, Q; kwargs...)

The Capture member whose polar motion sits on a double root of Θ(z), at constant latitude.
The constants must place z on that double root. Keywords as in
`kerr_geo_capture_component`.
"""
function kerr_geo_capture_constant_latitude(a::Real, energy::Real, lz::Real, q::Real;
        kwargs...)
    orbit = kerr_geo_capture_component(a, energy, lz, q;
        polar_sector=:constant_latitude, kwargs...)
    orbit.Status.polar.sector === :constant_latitude || error(
        "Constants do not select an exact constant-latitude polar double root.")
    return orbit
end

"""
    kerr_geo_capture_equator_attractive(a, E, Lz; kwargs...)

The Capture member with Q = 0 whose polar motion approaches the equator asymptotically
instead of lying in it (possible for E > 1). Keywords as in `kerr_geo_capture_component`.
"""
kerr_geo_capture_equator_attractive(a::Real, energy::Real, lz::Real; kwargs...) =
    kerr_geo_capture_component(a, energy, lz, 0.0; polar_sector=:equator_attractive,
        kwargs...)
