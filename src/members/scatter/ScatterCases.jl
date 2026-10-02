# Class D (Scatter) members D1, D2: component selection and assembly.

# Q = 0 scatter constants are taken as equatorial unless another sector is requested
_scatter_polar_sector(q, requested) =
    requested === nothing && iszero(q) ? :equatorial : requested

"""
    kerr_geo_scatter_component(a, E, Lz, Q; case_id=nothing, polar_sector=nothing,
                               polar_phase=0.0, polar_hemisphere=:north)

Classify (a, E, Lz, Q) and construct its Class D (Scatter) component: D1 (E = 1, three
real roots) or D2 (E > 1, four real roots), from infinity through the turning point (λ = 0,
where t, φ, τ are zero) back to infinity, with the closed-form r(λ) of `kerr_geo_scatter`. In
the pendular and equatorial sectors the polar motion is that of `kerr_geo_scatter`; in the
vortical, equator-attractive and axis-crossing ones that of the polar engine. Since
R(r₊) = P(r₊)² > 0 and R < 0 just inside the turning point, R has one more root between the
horizon and the turning point (r2 for D1, r3 for D2).
"""
function kerr_geo_scatter_component(a::Real, energy::Real, lz::Real, q::Real;
        case_id=nothing,
        polar_sector=nothing,
        polar_phase::Real=0.0,
        polar_hemisphere::Symbol=:north)
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=polar_sector)
    return _scatter_component(a, energy, lz, q, classification, case_id;
        polar_sector=polar_sector, polar_phase=polar_phase, polar_hemisphere=polar_hemisphere)
end

# the D1 or D2 member of classified constants
function _scatter_component(a, energy, lz, q, classification, case_id; polar_sector=nothing,
        polar_phase=0.0, polar_hemisphere::Symbol=:north)
    case_id === nothing || case_id in (:D1, :D2) || error(
        "`kerr_geo_scatter_component` builds D1 and D2, not $(case_id); D-H1 and D-H2 come from `kerr_geo_horizon_scatter`, D-X1 and D-X2 from `kerr_geo_extremal`.")
    classification = _with_polar_sector(classification, _scatter_polar_sector(q, polar_sector))
    component = _class_component(classification, :scatter, case_id)
    sector = component.PolarSector
    window = sector in (:equatorial, :pendular)
    window || sector in (:vortical, :equator_attractive, :axis_crossing) ||
        error("Class D polar motion is equatorial, pendular, vortical, equator-attractive or axis-crossing, not $(sector).")
    phase = float(polar_phase)
    polar = window ? _window_polar_motion(a, energy, lz, q, phase) :
        _polar_solution(a, energy, lz, q, sector, phase; hemisphere=polar_hemisphere)
    return _scatter_member(a, energy, lz, q, component.CaseId, polar,
        _root_radii(classification); component=component,
        formula_family=component.FormulaFamily)
end

# D1, D2 and the horizon-root D-H1, D-H2 (the horizon pole has zero residue): from infinity
# through the turning point (λ = 0, t = φ = τ = 0) back to infinity.
function _scatter_member(a, energy, lz, q, id, polar, radii; component=nothing,
        structure=component.Metadata.structure, formula_family, tier=kerr_geo_tier(id),
        status=(;), chart=nothing)
    formula = id in (:D1, :D_H1) ? :parabolic_scatter : :hyperbolic_scatter
    radial = _scatter_radial_model(formula, a, energy, lz, (roots=radii,))
    radial === nothing && error("The radial roots do not have the $(formula) structure.")
    λ_inf = radial.lambda_infinity
    check(λ) = (isfinite(λ) && -λ_inf < λ < λ_inf || throw(DomainError(λ,
        "Mino time must lie strictly between the infinity endpoints ±$(λ_inf).")); float(λ))
    radial_radius, radial_mino = radial.radius, radial.mino
    radius(λ) = radial_radius(abs(λ))
    # t, φ, τ: radial spectral engine (infinity → turning point → infinity) + polar primitive
    potential = _radial_potential_from_roots(a, energy, lz, q, structure)
    coords = _engine_coordinates(a, energy, lz, q, radius, _polar_primitive(polar);
        potential=potential, domain=(-λ_inf, λ_inf), ends=(:infinity, :infinity), turn=0.0, σ=1.0,
        λ_bl=0.0)
    function lambda_of_radius(r; branch=:outgoing)
        branch in (:incoming, :outgoing) || error("The branch must be :incoming or :outgoing.")
        λ = radial_mino(float(r))
        return branch === :incoming ? -λ : λ
    end
    track = (r=λ -> radius(check(λ)), check=check, check_bl=check, coords=coords,
        sign_r=λ -> sign(check(λ)),
        domain=(mino=(-λ_inf, λ_inf), endpoint_closed=(false, false),
            endpoint_roles=(:past_infinity, :future_infinity)),
        reference=(lambda0_event=:finite_turning_point, t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=radial.r_turn,
            tau_zero_event=:finite_turning_point, lambda_regular=nothing,
            polar_phase=polar.metadata.phase,
            polar_phase_convention=polar.metadata.phase_convention),
        # the increments are those of the outgoing leg
        trajectory=merge((lambda_of_radius=lambda_of_radius,),
            _radius_increments(coords, lambda_of_radius, 1.0; regular=false)))
    return _engine_member(:scatter, id, a, energy, lz, q, polar, track, component, structure,
        potential, tier, (radial=radial.radial,),
        merge((formula_family=formula_family, formula_kind=formula), status), chart)
end

function kerr_geo_scatter_component(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_scatter_component(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_scatter_component(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_scatter_component(a, constants...; kwargs...)
end

# the finite-window scatter orbit of the same constants (pendular and equatorial sectors)
function _scatter_window(kg::KerrGeoScatterComponent)
    kg.Component.PolarSector in (:equatorial, :pendular) || error(
        "Asymptotic diagnostics exist for the pendular and equatorial sectors only.")
    (; a, E, Lz, Q) = kg.ConstantsOfMotion
    return kerr_geo_scatter(a, E, Lz, Q; input=:constants,
        polar_phase=kg.ReferenceZero.polar_phase)
end

kerr_geo_scatter_asymptotic_diagnostics(kg::KerrGeoScatterComponent) =
    kerr_geo_scatter_asymptotic_diagnostics(_scatter_window(kg))
kerr_geo_scatter_asymptotic_state(kg::KerrGeoScatterComponent, branch::Symbol) =
    kerr_geo_scatter_asymptotic_state(_scatter_window(kg), branch)
