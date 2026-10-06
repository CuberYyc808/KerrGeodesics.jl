# Class B (Plunge) members B1-B9: analytic radial models and assembly.

# B1–B9: from the outer turning point `turn` (λ = 0, where t, φ, τ are zero) into the future
# horizon at λ_h (v = ψ = 0 there). `r_of(λ)` on [0, λ_h], `lambda_of_radius` its inverse.
function _plunge_member(a, energy, lz, q, component, polar, λ_h, turn, r_of, lambda_of_radius;
        roots, kind)
    λ_h > 0 || error("The horizon Mino time of $(component.CaseId) is not positive.")
    check, check_bl = _turning_to_horizon_checks(λ_h)
    potential = _radial_potential_from_roots(a, energy, lz, q, component.Metadata.structure)
    coords = _engine_coordinates(a, energy, lz, q, r_of, _polar_primitive(polar);
        potential=potential, domain=(0.0, λ_h), ends=(:turning, :horizon), σ=-1.0, λ_bl=0.0,
        λ_regular=λ_h, σ_regular=-1.0)
    track = (r=λ -> r_of(check(λ)), check=check, check_bl=check_bl, coords=coords,
        radial_track=r_of,
        sign_r=λ -> -1.0,
        domain=(mino=(0.0, λ_h), endpoint_closed=(true, true),
            endpoint_roles=(:finite_turning_point, :future_horizon), horizon_lambda=λ_h),
        reference=(lambda0_event=:finite_turning_point, t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=turn, tau_zero_event=:finite_turning_point,
            lambda_regular=λ_h, polar_phase=polar.metadata.phase,
            polar_phase_convention=polar.metadata.phase_convention),
        trajectory=merge((lambda_of_radius=lambda_of_radius,),
            _radius_increments(coords, lambda_of_radius, -1.0; regular=false)))
    return _engine_member(:plunge, component.CaseId, a, energy, lz, q, polar, track,
        component, component.Metadata.structure, potential, kerr_geo_tier(component.CaseId),
        (model=roots,), (formula_family=component.FormulaFamily, formula_kind=kind), nothing)
end

# Mino-time checks of a member from a turning point (λ = 0) into the horizon (λ = λ_h):
# (r, z, τ, v, ψ) on [0, λ_h], (t, φ) on [0, λ_h). λ within `_mino_endpoint_tol` below 0 is the
# turning point.
function _turning_to_horizon_checks(λ_h)
    function check(λ)
        -_mino_endpoint_tol(_float_type(λ, λ_h)) <= λ <= λ_h || throw(DomainError(λ, "Mino time must lie in [0, λ_h = $(λ_h)]."))
        return max(float(λ), 0.0)
    end
    check_bl(λ) = (check(λ) < λ_h || throw(DomainError(λ,
        "BL t and phi exclude the exact future-horizon endpoint.")); check(λ))
    return check, check_bl
end

# B1–B9: the closed-form radial model of the case (models/SimpleRootModels.jl,
# models/InteriorRepeatedRoots.jl), λ = 0 at the turning point. The polar convention of B1, B3
# and B4 is that of `kerr_geo_plunge` (z = 0 at λ = 0 moving north; `polar_phase` shifts the
# polar argument from there); the other cases take `polar_phase` as the phase at λ = 0.
function _plunge_radial_component(a, energy, lz, q, classification, component, polar_phase)
    id = component.CaseId
    structure = classification.Status.root_structure
    model = id in (:B7, :B8, :B9) ?
        interior_repeated_radial_model(id, a, energy, lz, q, structure) :
        _plunge_radial_model(id, energy, structure)
    phase = float(polar_phase)
    polar = if id in (:B1, :B3, :B4)
        roots = _polar_quadratic_roots(a, energy, lz, q)
        # k_θ = cu_small/cu_big and 1 − k_θ = √disc/cu_big
        shifted = _polar_solution(a, energy, lz, q, component.PolarSector,
            phase - _ellip_k(sqrt(max(roots.disc, 0.0)) / roots.cu_big))
        merge(shifted, (metadata=merge(shifted.metadata, (phase=phase,
            phase_convention=:equator_northward_at_zero_phase)),))
    else
        _polar_solution(a, energy, lz, q, component.PolarSector, phase)
    end
    rplus = kerr_horizons(a).rplus
    relative = _horizon_relative_model(a,energy,lz,q,model)
    λ_h = relative.horizon_time
    turn = model.turn
    r_of = relative
    function lambda_of_radius(r)
        rplus <= r <= turn || throw(DomainError(r,
            "The radius lies outside [r+, r_turn] of $(id)."))
        r==rplus && return λ_h
        r==turn && return 0.0
        gap=_wide_sub(_wide(float(r)),relative.horizon)
        return _horizon_lambda_of_gap(relative,gap[1]+gap[2])
    end
    return _plunge_member(a, energy, lz, q, component, polar, λ_h, model.turn, r_of,
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
    return _plunge_component(a, energy, lz, q, classification, case_id;
        polar_phase=polar_phase, reference_radius=reference_radius)
end

# the Plunge member `case_id` of classified constants
function _plunge_component(a, energy, lz, q, classification, case_id; polar_phase=0.0,
        reference_radius=nothing)
    component = _class_component(classification, :plunge, case_id)
    # the turning point is the reference event of every Plunge member
    reference_radius === nothing || isapprox(reference_radius,
        component.UpperEndpoint.Radius;
        atol=_radius_tol(typeof(component.UpperEndpoint.Radius)),
        rtol=_radius_tol(typeof(component.UpperEndpoint.Radius))) || error(
        "Plunge members use their turning point as the reference radius.")
    return _plunge_radial_component(a, energy, lz, q, classification, component, polar_phase)
end

function kerr_geo_plunge_component(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_plunge_component(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_plunge_component(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_plunge_component(a, constants...; kwargs...)
end
