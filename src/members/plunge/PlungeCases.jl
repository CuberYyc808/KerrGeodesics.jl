# Class B (Plunge) members B1-B9: analytic radial models and assembly.

const _RADIUS_TOL = 2.0e-11

function _b2_model(energy, roots)
    x1, x2, x3 = roots
    kappa = -_e2m1(energy)
    b = x2 - x1
    c = x3 - x1
    n = b / c
    scale = 2 / (sqrt(kappa) * c)
    function phi_of_r(r)
        return asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))
    end
    function basis(r)
        phi = phi_of_r(r)
        pin = _pi_real(n, phi, 0.0)
        i0 = scale * pin
        i1 = scale * (x3 * pin - c * phi)
        i2 = scale * (x3^2 * pin - 2 * x3 * c * phi +
            c^2 * ((1 - n / 2) * phi + n * sin(2 * phi) / 4))
        return (I0=i0, I1=i1, I2=i2)
    end
    function pole(h, r)
        phi = phi_of_r(r)
        nh = b / (h - x1)
        return scale / (x1 - h) * (
            -n / (nh - n) * _pi_real(n, phi, 0.0) +
            nh / (nh - n) * _pi_real(nh, phi, 0.0))
    end
    function inverse_i0(target)
        k = sqrt(1 - n)
        angle = clamp(k * target / scale, 0.0, pi / 2)
        angle >= pi / 2 - 8 * eps(Float64) && return x2
        tangent = tan(angle) / k
        s2 = tangent^2 / (1 + tangent^2)
        return x1 + b * s2
    end
    return (kind=:b2_outer_double_inner_simple, roots=(x1=x1, x2=x2, x3=x3),
        upper=x2, basis=basis, pole=pole,
        inverse_i0=inverse_i0)
end

function _b5_model(roots)
    x1, x2, x3 = roots
    b = x2 - x1
    c = x3 - x1
    m = b / c
    scale = sqrt(2 / c)
    phi_of_r(r) = asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))
    function j2(phi)
        return (Elliptic.F(phi, m) - Elliptic.E(phi, m)) / m
    end
    function j4(phi)
        boundary = sin(phi) * cos(phi) *
            sqrt(max(1 - m * sin(phi)^2, 0.0))
        return (boundary - Elliptic.F(phi, m) +
            (2 + 2 * m) * j2(phi)) / (3 * m)
    end
    function basis(r)
        phi = phi_of_r(r)
        f = Elliptic.F(phi, m)
        return (
            I0=scale * f,
            I1=scale * (x1 * f + b * j2(phi)),
            I2=scale * (x1^2 * f + 2 * x1 * b * j2(phi) + b^2 * j4(phi)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        return scale * _pi_real(b / (h - x1), phi, m) / (x1 - h)
    end
    function inverse_i0(target)
        sn = Elliptic.Jacobi.sn(clamp(target / scale, 0.0, Elliptic.K(m)), m)
        return x1 + b * sn^2
    end
    return (kind=:b5_parabolic_three_simple, roots=(x1=x1, x2=x2, x3=x3),
        upper=x2, basis=basis, pole=pole,
        inverse_i0=inverse_i0)
end

function _b6_model(energy, roots)
    x1, x2, x3, x4 = roots
    lead = _e2m1(energy)
    n = (x3 - x2) / (x4 - x2)
    m = (x3 - x2) * (x4 - x1) / ((x4 - x2) * (x3 - x1))
    g = x3 - x4
    scale = 2 / sqrt(lead * (x4 - x2) * (x3 - x1))
    function phi_of_r(r)
        s2 = (x4 - x2) * (x3 - r) / ((x3 - x2) * (x4 - r))
        return asin(sqrt(clamp(s2, 0.0, 1.0)))
    end
    function basis(r)
        phi = phi_of_r(r)
        f = Elliptic.F(phi, m)
        pin = _pi_real(n, phi, m)
        return (
            I0=-scale * f,
            I1=-scale * (x4 * f + g * pin),
            I2=-scale * (x4^2 * f + 2 * x4 * g * pin +
                g^2 * _j2_legendre(n, m, phi)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        nh = n * (x4 - h) / (x3 - h)
        return -scale * (Elliptic.F(phi, m) / (x4 - h) +
            (x4 - x3) * _pi_real(nh, phi, m) / ((x4 - h) * (x3 - h)))
    end
    function inverse_i0(target)
        f = clamp(-target / scale, 0.0, Elliptic.K(m))
        s2 = Elliptic.Jacobi.sn(f, m)^2
        numerator = (x4 - x2) * x3 - s2 * (x3 - x2) * x4
        denominator = (x4 - x2) - s2 * (x3 - x2)
        return numerator / denominator
    end
    return (kind=:b6_hyperbolic_four_simple, roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        upper=x3, basis=basis, pole=pole,
        inverse_i0=inverse_i0)
end

function _plunge_analytic_model(case_id, energy, radii)
    case_id === :B2 && return _b2_model(energy, (radii[1], radii[2], radii[3]))
    case_id === :B5 && return _b5_model((radii[1], radii[2], radii[3]))
    case_id === :B6 && return _b6_model(energy, (radii[1], radii[2], radii[3], radii[4]))
    error("The closed-form Plunge models here are those of B2, B5 and B6, not $(case_id).")
end

# B1–B9: from the outer turning point `turn` (λ = 0, where t, φ, τ are zero) into the future
# horizon at λ_h (v = ψ = 0 there). `r_of(λ)` on [0, λ_h], `lambda_of_radius` its inverse.
function _plunge_member(a, energy, lz, q, component, polar, λ_h, turn, r_of, lambda_of_radius;
        roots, kind)
    λ_h > 0 || error("The horizon Mino time of $(component.CaseId) is not positive.")
    check, check_bl = _turning_to_horizon_checks(λ_h)
    coords = _engine_coordinates(a, energy, lz, q, r_of, _polar_primitive(polar);
        domain=(0.0, λ_h), ends=(:turning, :horizon), σ=-1.0, λ_bl=0.0, λ_regular=λ_h,
        σ_regular=-1.0)
    track = (r=λ -> r_of(check(λ)), check=check, check_bl=check_bl, coords=coords,
        sign_r=λ -> -1.0,
        domain=(mino=(0.0, λ_h), endpoint_closed=(true, true),
            endpoint_roles=(:finite_turning_point, :future_horizon), horizon_lambda=λ_h),
        reference=(lambda0_event=:finite_turning_point, t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=turn, tau_zero_event=:finite_turning_point,
            lambda_regular=λ_h, polar_phase=polar.metadata.phase,
            polar_phase_convention=polar.metadata.phase_convention),
        trajectory=merge((lambda_of_radius=lambda_of_radius,),
            _radius_increments(coords, lambda_of_radius, -1.0; regular=false)))
    return _engine_member(:plunge, component.CaseId, a, energy, lz, q, polar, track;
        component=component, roots=(model=roots,),
        status=(formula_family=component.FormulaFamily, formula_kind=kind))
end

# Mino-time checks of a member from a turning point (λ = 0) into the horizon (λ = λ_h):
# (r, z, τ, v, ψ) on [0, λ_h], (t, φ) on [0, λ_h). λ within 2e-12 below 0 is the turning point.
function _turning_to_horizon_checks(λ_h)
    function check(λ)
        -2.0e-12 <= λ <= λ_h || throw(DomainError(λ, "Mino time must lie in [0, λ_h = $(λ_h)]."))
        return max(float(λ), 0.0)
    end
    check_bl(λ) = (check(λ) < λ_h || throw(DomainError(λ,
        "BL t and phi exclude the exact future-horizon endpoint.")); check(λ))
    return check, check_bl
end

"""
B1, B3, B4: r(λ) from the shared closed-form models (models/SimpleRootModels.jl). Their polar
convention is the one of `kerr_geo_plunge`: z = 0 at λ = 0 moving north (`polar_phase` shifts
the polar argument from there).
"""
function _plunge_model_component(a, energy, lz, q, component, polar_phase)
    model = _plunge_radial_model(component.CaseId, a, energy, lz, q)
    rplus = kerr_horizons(a).rplus
    _, _, k_theta = _plunge_polar_parameters(a, energy, lz, q)
    polar = _polar_solution(a, energy, lz, q, component.PolarSector,
        float(polar_phase) - Elliptic.K(k_theta))
    polar = merge(polar, (metadata=merge(polar.metadata, (phase=float(polar_phase),
        phase_convention=:equator_northward_at_zero_phase)),))
    λ_h = model.lambda_horizon
    r_of(λ) = λ >= λ_h ? rplus : clamp(model.r(λ), rplus, model.turn)
    function lambda_of_radius(r)
        rplus <= r <= model.turn || throw(DomainError(r,
            "The radius lies outside [r+, r_turn] of $(component.CaseId)."))
        return r == rplus ? λ_h : model.lambda_of_r(r)
    end
    return _plunge_member(a, energy, lz, q, component, polar, λ_h, model.turn, r_of,
        lambda_of_radius; roots=model.roots, kind=model.kind)
end

# B2, B5, B6 (closed-form models) and B7–B9 (repeated roots below the horizon).
function _analytic_component(a, energy, lz, q, classification, component, polar_phase)
    rplus = kerr_horizons(a).rplus
    model = component.CaseId in (:B7, :B8, :B9) ?
        interior_repeated_radial_model(component.CaseId, a, energy, lz, q,
            classification.Status.root_structure) :
        _plunge_analytic_model(component.CaseId, energy, _root_radii(classification))
    polar = _polar_solution(a, energy, lz, q, component.PolarSector, polar_phase)
    c = model.basis(model.upper).I0
    λ_h = c - model.basis(rplus).I0
    r_of(λ) = λ >= λ_h ? rplus : clamp(model.inverse_i0(c - λ), rplus, model.upper)
    function lambda_of_radius(r)
        rplus <= r <= model.upper || throw(DomainError(r,
            "The radius lies outside [r+, r_turn] of $(component.CaseId)."))
        return r == rplus ? λ_h : c - model.basis(r).I0
    end
    return _plunge_member(a, energy, lz, q, component, polar, λ_h, model.upper, r_of,
        lambda_of_radius; roots=model.roots, kind=model.kind)
end

"""
    kerr_geo_plunge_component(a, E, Lz, Q; case_id=nothing, polar_sector=nothing,
                              polar_phase=0.0, reference_radius=nothing)

The Plunge member `case_id` (B1–B9) of the constants `(a, E, Lz, Q)`, a
`KerrGeoPlungeComponent`: from its outer turning point into the future horizon. λ = 0 at the
turning point, where t, φ and τ vanish; v = ψ = 0 on the future horizon. `polar_phase` is the
polar phase at λ = 0 (convention in `ReferenceZero.polar_phase_convention`);
`reference_radius`, when given, must be the turning point. r(λ) comes from the closed-form
radial model of each case, and t, φ, τ, v, ψ from the spectral coordinate engine. The plunges
from a Critical root (K2, K5, K8, K11) are built by `kerr_geo_critical_component`.
"""
function kerr_geo_plunge_component(a::Real, energy::Real, lz::Real, q::Real;
        case_id=nothing,
        polar_sector=nothing,
        polar_phase::Real=0.0,
        reference_radius=nothing)
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=polar_sector)
    component = _class_component(classification, :plunge, case_id)
    # the turning point is the reference event of every Plunge member
    reference_radius === nothing || isapprox(reference_radius,
        component.UpperEndpoint.Radius; atol=_RADIUS_TOL, rtol=_RADIUS_TOL) || error(
        "Plunge members use their turning point as the reference radius.")
    component.CaseId in (:B1, :B3, :B4) &&
        return _plunge_model_component(a, energy, lz, q, component, polar_phase)
    return _analytic_component(a, energy, lz, q, classification, component,
        float(polar_phase))
end

function kerr_geo_plunge_component(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_plunge_component(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_plunge_component(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_plunge_component(a, constants...; kwargs...)
end
