# Class C (Capture) members C1–C12: from infinity into the future horizon. Radial models of C2
# and C4, the closed-form C1/C3 member, the shared model-based member (C2, C4, C5, C6–C12) and
# the exported constructors of the polar sectors.

const _NONGENERIC_POLAR_SECTORS = (
    :vortical, :constant_latitude, :equator_attractive, :axis_crossing)

function _c2_model(roots)
    x1, x2, x3 = roots
    aroot = x3 - x1
    m = (x2 - x1) / aroot
    scale = sqrt(2 / aroot)
    function phi_of_r(r)
        return asin(sqrt(clamp(aroot / (r - x1), 0.0, 1.0)))
    end
    w(phi) = sqrt(max(1 - m * sin(phi)^2, 0.0))
    function j2(phi)
        return -cos(phi) / sin(phi) * w(phi) +
            Elliptic.F(phi, m) - Elliptic.E(phi, m)
    end
    function j4(phi)
        return ((2 + 2 * m) * j2(phi) - m * Elliptic.F(phi, m) -
            cos(phi) / sin(phi) * w(phi) / sin(phi)^2) / 3
    end
    function basis(r)
        phi = phi_of_r(r)
        f = Elliptic.F(phi, m)
        return (
            I0=-scale * f,
            I1=-scale * (x1 * f + aroot * j2(phi)),
            I2=-scale * (x1^2 * f + 2 * x1 * aroot * j2(phi) +
                aroot^2 * j4(phi)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        n = (h - x1) / aroot
        abs(n) > 1.0e-14 || error(
            "The C2 pole primitive is undefined at zero characteristic (h = x1), where the horizon residue vanishes.")
        return -scale * (_pi_real(n, phi, m) - Elliptic.F(phi, m)) /
            (aroot * n)
    end
    function inverse_i0(target)
        f = max(-target / scale, 0.0)
        f > 0 || return Inf
        sn = Elliptic.Jacobi.sn(f, m)
        abs(sn) > 0 || return Inf
        return x1 + aroot / sn^2
    end
    return (
        kind=:c2_parabolic_three_real_direct_capture,
        roots=(x1=x1, x2=x2, x3=x3),
        lower=x3,
        basis=basis,
        pole=pole,
        inverse_i0=inverse_i0,
        infinity_i0=0.0,
    )
end

function _c4_model(energy, roots)
    x1, x2, x3, x4 = roots
    lead = _e2m1(energy)
    m = (x4 - x1) * (x3 - x2) / ((x4 - x2) * (x3 - x1))
    n = (x4 - x1) / (x3 - x1)
    scale = 2 / sqrt(lead * (x4 - x2) * (x3 - x1))
    function phi_of_r(r)
        s2 = (x3 - x1) * (r - x4) / ((x4 - x1) * (r - x3))
        return asin(sqrt(clamp(s2, 0.0, 1.0)))
    end
    function basis(r)
        phi = phi_of_r(r)
        f = Elliptic.F(phi, m)
        pin = _pi_real(n, phi, m)
        return (
            I0=scale * f,
            I1=scale * (x3 * f + (x4 - x3) * pin),
            I2=scale * (x3^2 * f + 2 * x3 * (x4 - x3) * pin +
                (x4 - x3)^2 * _j2_legendre(n, m, phi)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        nh = (x4 - x1) * (x3 - h) / ((x3 - x1) * (x4 - h))
        return scale * (Elliptic.F(phi, m) / (x3 - h) +
            (x3 - x4) * _pi_real(nh, phi, m) /
            ((x4 - h) * (x3 - h)))
    end
    function inverse_i0(target)
        f = clamp(target / scale, 0.0, Elliptic.K(m))
        s2 = Elliptic.Jacobi.sn(f, m)^2
        alpha = x3 - x1
        beta = x4 - x1
        denominator = alpha - s2 * beta
        abs(denominator) > 1.0e-15 || return Inf
        return (alpha * x4 - s2 * beta * x3) / denominator
    end
    phi_infinity = asin(sqrt((x3 - x1) / (x4 - x1)))
    return (
        kind=:c4_hyperbolic_four_real_direct_capture,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        lower=x4,
        basis=basis,
        pole=pole,
        inverse_i0=inverse_i0,
        infinity_i0=scale * Elliptic.F(phi_infinity, m),
    )
end

function _capture_analytic_model(case_id, energy, radii)
    case_id === :C2 && return _c2_model((radii[1], radii[2], radii[3]))
    case_id === :C4 && return _c4_model(
        energy, (radii[1], radii[2], radii[3], radii[4]))
    error("The analytic radial model covers C2 and C4, not $(case_id).")
end

# C1 (E = 1, one real root) and C3 (E > 1, two real roots): the closed-form radial formulas
# of the finite-window capture API. λ = 0 on the future horizon, where τ = v = ψ = 0; t and φ
# are zero at the reference radius. The pendular and equatorial sectors take the polar motion
# of `kerr_geo_capture`, the other sectors that of the polar engine.
function _closed_form_capture_component(a, energy, lz, q, component, polar_phase,
        polar_hemisphere, reference_radius)
    if component.CaseId === :C1
        radial = _c1_one_real_parameters(a, lz, q)
        radial === nothing && error(
            "The radial roots do not have the C1 structure (one real root below r₊ and a complex pair).")
        return _closed_form_capture(a, energy, lz, q, component, polar_phase,
            polar_hemisphere, reference_radius, radial.rplus,
            r -> _c1_lambda_of_r(radial, r), _c1_lambda_infinity(radial),
            λ -> _c1_radius_from_horizon_lambda(radial, λ),
            (real=radial.real_roots, complex=(rho=radial.rho, eta=radial.eta)),
            :c1_parabolic_one_real_complex)
    end
    radial = _c3_complex_parameters(a, energy, lz, q)
    radial === nothing && error(
        "The radial roots do not have the C3 structure (two real roots below r₊ and a complex pair).")
    return _closed_form_capture(a, energy, lz, q, component, polar_phase, polar_hemisphere,
        reference_radius, radial.rplus, r -> _c3_lambda_of_r(radial, r),
        _c3_lambda_infinity(radial), λ -> _c3_radius_from_horizon_lambda(radial, λ),
        (real=radial.real_roots, complex=(rho=radial.rho, eta=radial.eta)),
        :c3_hyperbolic_two_real_complex)
end

function _closed_form_capture(a, energy, lz, q, component, polar_phase, polar_hemisphere,
        reference_radius, rplus, λ_of_r, λ_span, r_of, radial_roots, formula_kind)
    λ_h = λ_of_r(rplus)
    λ_inf = λ_h - λ_span
    function lambda_of_radius(r)
        r >= rplus || throw(DomainError(r, "The radius must lie outside the future horizon."))
        return λ_h - λ_of_r(r)
    end
    check, check_bl = _infinity_to_horizon_checks(λ_inf)
    window = !(component.PolarSector in _NONGENERIC_POLAR_SECTORS)
    polar = window ? _window_polar_motion(a, energy, lz, q, polar_phase) :
        _polar_solution(a, energy, lz, q, component.PolarSector, polar_phase;
            hemisphere=polar_hemisphere)
    λ_ref = reference_radius === nothing ? 0.5 * λ_inf : lambda_of_radius(float(reference_radius))
    λ_inf < λ_ref < 0 || error("The reference radius must be a finite exterior point.")
    coords = _engine_coordinates(a, energy, lz, q, r_of, _polar_primitive(polar);
        domain=(λ_inf, 0.0), ends=(:infinity, :horizon), σ=-1.0, λ_bl=λ_ref,
        λ_regular=0.0, σ_regular=-1.0)
    track = (r=λ -> r_of(check(λ)), check=check, check_bl=check_bl, coords=coords,
        sign_r=λ -> -1.0,
        domain=(mino=(λ_inf, 0.0), endpoint_closed=(false, true),
            endpoint_roles=(:past_infinity, :future_horizon), horizon_lambda=0.0),
        reference=(lambda0_event=:future_horizon, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=λ_ref, t_phi_zero_radius=r_of(λ_ref),
            tau_zero_event=:future_horizon, lambda_regular=0.0,
            polar_phase=polar.metadata.phase,
            polar_phase_convention=polar.metadata.phase_convention),
        trajectory=merge((lambda_of_radius=lambda_of_radius,),
            _radius_increments(coords, lambda_of_radius, -1.0)))
    return _engine_member(:capture, component.CaseId, a, energy, lz, q, polar, track;
        component=component, roots=(model=radial_roots,),
        status=(formula_family=component.FormulaFamily, formula_kind=formula_kind))
end

# C2, C4, C5, C6–C12: a radial model from infinity into the future horizon (λ = 0, where
# τ = v = ψ = 0); t and φ are zero at the reference radius (default: λ = λ_∞/2). All radial
# increments (t, φ, τ, v, ψ) come from the spectral engine.
function _direct_component(a, energy, lz, q, component, model, polar, reference_radius;
        roots=model.roots, status=(;))
    rplus = kerr_horizons(a).rplus
    c = model.basis(rplus).I0
    λ_inf = c - model.infinity_i0
    λ_inf < 0 || error("The infinity endpoint must precede the future horizon.")
    check, check_bl = _infinity_to_horizon_checks(λ_inf)
    r_of(λ) = λ == 0 ? rplus : model.inverse_i0(c - λ)
    function lambda_of_radius(r)
        r > rplus || throw(DomainError(r, "The radius must lie outside the future horizon."))
        return c - model.basis(r).I0
    end
    λ_ref = reference_radius === nothing ? 0.5 * λ_inf : lambda_of_radius(float(reference_radius))
    λ_inf < λ_ref < 0 || error("The reference radius must be a finite exterior point.")
    coords = _engine_coordinates(a, energy, lz, q, r_of, _polar_primitive(polar);
        domain=(λ_inf, 0.0), ends=(:infinity, :horizon), σ=-1.0, λ_bl=λ_ref,
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
    return _engine_member(:capture, component.CaseId, a, energy, lz, q, polar, track;
        component=component, roots=(model=roots,),
        status=merge((formula_family=component.FormulaFamily, formula_kind=model.kind),
            status))
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
    component = _class_component(classification, :capture, case_id)
    phase = float(polar_phase)
    id = component.CaseId
    id === :C5 && return kerr_geo_capture_four_complex(a, energy, lz, q;
        polar_sector=component.PolarSector, polar_hemisphere=polar_hemisphere,
        polar_phase=phase, reference_radius=reference_radius)
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
