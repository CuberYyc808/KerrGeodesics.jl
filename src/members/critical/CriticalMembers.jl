# Class K (Critical) members: radial motion on, or tending asymptotically to, an unstable or
# marginal repeated root r_c > r₊ (see `kerr_geo_is_critical`). The members of one set of
# constants lie on the root (K1, K3, K6, K9), on its outer side (K4 homoclinic, K7/K10 from
# infinity) or on its inner side into the future horizon (K2, K5, K8, K11); the triple root K1
# pairs with K2 only. One assembly builds them all; a role only supplies its radial track
# (r(λ), λ(r) and the engine endpoints) and the radial models are in Models.jl.

"""
    CRITICAL_ROLES

The Critical cases by role, the side of the repeated root on which the member lies:
`on_root = (:K1, :K3, :K6, :K9)`, `outer = (:K4, :K7, :K10)` and
`inner = (:K2, :K5, :K8, :K11)`. The role is a position, not a direction of motion: K7 and
K10 move inward on the outer side.
"""
const CRITICAL_ROLES = (on_root=(:K1, :K3, :K6, :K9), outer=(:K4, :K7, :K10),
    inner=(:K2, :K5, :K8, :K11))

"""
    kerr_geo_critical_role(case_id)

The role of a Critical case, the side of the repeated root it lies on: `:on_root`, `:outer`
or `:inner` (see `CRITICAL_ROLES`).
"""
function kerr_geo_critical_role(case_id::Symbol)
    for role in keys(CRITICAL_ROLES)
        case_id in CRITICAL_ROLES[role] && return role
    end
    throw(ArgumentError("$(case_id) is not a Critical case."))
end

"""
    kerr_geo_critical_component(a, E, Lz, Q; case_id, polar_sector=nothing, polar_phase=0.0,
                                polar_hemisphere=:north, reference_radius=nothing)

The Critical member `case_id` (K1–K11) of the constants `(a, E, Lz, Q)`, a
`KerrGeoCriticalComponent`; `rc` below is the repeated root.

Zero points (`ReferenceZero`):
- on the root (K1, K3, K6, K9) and on the homoclinic orbit K4 (at apastron), t, φ and τ
  vanish at λ = 0;
- from infinity (K7, K10), λ = 0 and t = φ = τ = v = ψ = 0 at `reference_radius`
  (default `rc + max(1, rc − r₊)`);
- into the horizon (K2, K5, K8, K11), λ = 0, τ = 0 and v = ψ = 0 on the future horizon, and
  t = φ = 0 at `reference_radius` (default `r₊ + 0.55(rc − r₊)`, `(r₊ + rc)/2` for K5).

`polar_phase` is the polar phase at λ = 0, except for K2, K8 and K11, where it is the phase
at the reference-radius event (`ReferenceZero.polar_phase_event = :reference_radius`).
"""
function kerr_geo_critical_component(a::Real, energy::Real, lz::Real, q::Real;
        case_id::Symbol, polar_sector=nothing, kwargs...)
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=polar_sector)
    return _critical_component(a, energy, lz, q, classification, case_id;
        polar_sector=polar_sector, kwargs...)
end

kerr_geo_critical_component(a::Real, constants::NamedTuple; kwargs...) =
    kerr_geo_critical_component(a, constants.E, constants.Lz, constants.Q; kwargs...)
kerr_geo_critical_component(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...) =
    kerr_geo_critical_component(a, constants...; kwargs...)

# the member `case_id` of classified constants
function _critical_component(a, energy, lz, q, classification, case_id;
        polar_sector=nothing, polar_phase=0.0, polar_hemisphere::Symbol=:north,
        reference_radius=nothing)
    component = _class_component(classification, :critical, case_id)
    a, energy, lz, q = float(a), float(energy), float(lz), float(q)
    phase = float(polar_phase)
    sector = _resolved_polar_sector(classification, polar_sector, energy, lz, q)
    polar = _polar_solution(a, energy, lz, q, sector, phase; hemisphere=polar_hemisphere)
    root = component.LowerEndpoint.Multiplicity >= 2 ? component.LowerEndpoint :
        component.UpperEndpoint
    rc = float(root.Radius)
    track = case_id in CRITICAL_ROLES.on_root ? _on_root_track(a, energy, lz, polar, rc) :
        _radial_track(a, energy, lz, q, case_id, classification, polar, rc, reference_radius)
    return _critical_member(a, energy, lz, q, case_id, component, classification,
        get(track, :polar, polar), phase, root, track)
end

# the member of one role of these constants
function _critical_role_member(role, a, energy, lz, q; polar_sector=nothing, kwargs...)
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=polar_sector)
    ids = [c.CaseId for c in classification.Components if c.CaseId in CRITICAL_ROLES[role]]
    isempty(ids) && error("These constants have no Critical member of role $(role); " *
        "their cases are $(classification.CaseIds).")
    return _critical_component(a, energy, lz, q, classification, only(ids);
        polar_sector=polar_sector, kwargs...)
end

"""
    kerr_geo_critical_spherical(a, E, Lz, Q; kwargs...)

The Critical member that sits on the repeated root: the ISCO or ISSO (K1), or an unstable
circular or spherical orbit with E < 1 (K3), E = 1 (K6) or E > 1 (K9). Keywords as in
`kerr_geo_critical_component`.
"""
kerr_geo_critical_spherical(a::Real, energy::Real, lz::Real, q::Real; kwargs...) =
    _critical_role_member(:on_root, a, energy, lz, q; kwargs...)

"""
    kerr_geo_critical_homoclinic(a, E, Lz, Q; kwargs...)

The homoclinic orbit K4 (E < 1): it leaves the unstable spherical orbit, turns at apastron
(λ = 0) and returns to the same orbit. Keywords as in `kerr_geo_critical_component`.
"""
kerr_geo_critical_homoclinic(a::Real, energy::Real, lz::Real, q::Real; kwargs...) =
    kerr_geo_critical_component(a, energy, lz, q; case_id=:K4, kwargs...)

"""
    kerr_geo_critical_plunge(a, E, Lz, Q; kwargs...)

The Critical member that leaves the repeated root inward and crosses the future horizon:
K2, K5, K8 or K11. Keywords as in `kerr_geo_critical_component`.
"""
kerr_geo_critical_plunge(a::Real, energy::Real, lz::Real, q::Real; kwargs...) =
    _critical_role_member(:inner, a, energy, lz, q; kwargs...)

# ---- radial tracks (see Assembly.jl for the fields of a track) --------------------------------

# On the root, t, φ, τ grow at the constant radial rates plus the polar primitives.
function _on_root_track(a, energy, lz, polar, rc)
    momentum = kerr_radial_momentum(a, energy, lz, rc)
    rate_t = (rc^2 + a^2) * momentum / kerr_delta(a, rc)
    rate_phi = a * momentum / kerr_delta(a, rc) - a * energy
    p = _polar_primitive(polar)
    coords = (t=λ -> rate_t * λ + p(λ)[1], phi=λ -> rate_phi * λ + p(λ)[2],
        tau=λ -> rc^2 * λ + p(λ)[3])
    return (r=λ -> (_finite_mino(λ); rc), check=_finite_mino, check_bl=_finite_mino,
        coords=coords, sign_r=λ -> 0.0, model=nothing,
        domain=(mino=(-Inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:infinite_past_worldline, :infinite_future_worldline)),
        reference=(lambda0_event=:polar_phase_reference, t_phi_zero_event=:polar_phase_reference,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=rc, tau_zero_event=:polar_phase_reference,
            lambda_regular=nothing),
        trajectory=(;))
end

function _radial_track(a, energy, lz, q, id, classification, polar, rc, reference_radius)
    radii = id in (:K4, :K5) ? Tuple(_double_root_factorization(a, energy, lz, q, rc)) :
        _root_radii(classification)
    model = _critical_radial_model(id, energy, radii)
    rplus = kerr_horizons(a).rplus
    track = id === :K4 ? _homoclinic_track(a, energy, lz, q, model, polar) :
        id in CRITICAL_ROLES.outer ?
        _infinity_track(a, energy, lz, q, model, polar, rplus, reference_radius) :
        _inward_track(a, energy, lz, q, model, polar, rplus, reference_radius, id)
    # radial increments between radii of one monotone leg, from the engine: the outgoing leg of
    # K4 (no regular chart), the incoming ones of the others (with v, ψ)
    λ_of = id === :K4 ? (r -> track.trajectory.lambda_of_radius(r; branch=:outgoing)) :
        track.trajectory.lambda_of_radius
    increments = _radius_increments(track.coords, λ_of, id === :K4 ? 1.0 : -1.0;
        regular=id !== :K4)
    return merge(track, (model=model, trajectory=merge(track.trajectory, increments)))
end

# K4: λ = 0 at the apastron r_a and r(λ) is even in λ; both ends tend to r_c.
function _homoclinic_track(a, energy, lz, q, model, polar)
    radius(λ) = model.inverse_i0(-abs(_finite_mino(λ)))
    coords = _engine_coordinates(a, energy, lz, q, radius, _polar_primitive(polar);
        domain=(-Inf, Inf), ends=(:asymptote, :asymptote), turn=0.0, σ=-1.0,
        rd=model.lower, λ_bl=0.0)
    return (r=radius, check=_finite_mino, check_bl=_finite_mino, coords=coords,
        sign_r=λ -> -sign(_finite_mino(λ)),
        domain=(mino=(-Inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:past_repeated_root_asymptote, :future_repeated_root_asymptote)),
        reference=(lambda0_event=:finite_turning_point, t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=model.upper,
            tau_zero_event=:finite_turning_point, lambda_regular=nothing),
        # λ < 0 is the outgoing leg (r_c → r_a), λ > 0 the incoming one
        trajectory=(lambda_of_radius=(r; branch=:outgoing) -> (i = model.basis(r).I0;
            branch === :outgoing ? i : branch === :incoming ? -i :
            error("The branch must be :incoming or :outgoing.")),))
end

# K7, K10: from infinity towards r_c; λ = 0 (and t, φ, τ, v, ψ = 0) at the reference radius.
function _infinity_track(a, energy, lz, q, model, polar, rplus, reference_radius)
    rc = model.lower
    rref = reference_radius === nothing ? rc + max(1.0, rc - rplus) : float(reference_radius)
    rref > rc || error("The reference radius must lie above the repeated root.")
    c = model.basis(rref).I0
    λ_inf = c - model.infinity_i0
    check(λ) = (_finite_mino(λ) > λ_inf || throw(DomainError(λ,
        "Mino time must follow the infinity endpoint λ = $(λ_inf).")); float(λ))
    radius(λ) = model.inverse_i0(c - check(λ))
    function lambda_of_radius(r)
        r > rc || throw(DomainError(r, "The radius must lie above the repeated root."))
        return c - model.basis(r).I0
    end
    coords = _engine_coordinates(a, energy, lz, q, radius, _polar_primitive(polar);
        domain=(λ_inf, Inf), ends=(:infinity, :asymptote), σ=-1.0, rd=rc,
        λ_bl=0.0, λ_regular=0.0, σ_regular=-1.0)
    return (r=radius, check=check, check_bl=check, coords=coords, sign_r=λ -> -1.0,
        domain=(mino=(λ_inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:past_infinity, :future_repeated_root_asymptote)),
        reference=(lambda0_event=:reference_radius, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=rref, tau_zero_event=:reference_radius,
            lambda_regular=0.0),
        trajectory=(lambda_of_radius=lambda_of_radius,))
end

# K2, K5, K8, K11: from r_c into the future horizon, λ = 0 on the horizon (where v = ψ = 0 and
# τ = 0); t, φ vanish at the reference radius. The polar phase of K2, K8, K11 is given at the
# reference-radius event (the polar solution is delayed to it); K5's phase is given on the horizon.
function _inward_track(a, energy, lz, q, model, polar, rplus, reference_radius, id)
    rc = model.upper
    rref = reference_radius !== nothing ? float(reference_radius) :
        id === :K5 ? 0.5 * (rplus + rc) : rplus + 0.55 * (rc - rplus)
    rplus < rref < rc || error("The reference radius lies outside (r+, r_c).")
    c = model.basis(rplus).I0
    λ_h = 0.0
    check(λ) = (_finite_mino(λ) <= λ_h || throw(DomainError(λ,
        "Mino time lies beyond the future-horizon endpoint λ = $(λ_h).")); float(λ))
    check_bl(λ) = (check(λ) < λ_h || throw(DomainError(λ,
        "BL t and phi exclude the exact future-horizon endpoint.")); float(λ))
    radius(λ) = (λ = check(λ); λ == λ_h ? rplus : clamp(model.inverse_i0(c - λ), rplus, rc))
    function lambda_of_radius(r)
        rplus <= r <= rc || throw(DomainError(r, "The radius lies outside [r+, r_c]."))
        return r == rplus ? λ_h : c - model.basis(r).I0
    end
    λ_ref = c - model.basis(rref).I0
    polar = id === :K5 ? polar : _polar_delayed(polar, λ_ref)
    coords = _engine_coordinates(a, energy, lz, q, radius, _polar_primitive(polar);
        domain=(-Inf, λ_h), ends=(:asymptote, :horizon), σ=-1.0, rd=rc,
        λ_bl=λ_ref, λ_regular=λ_h, σ_regular=-1.0)
    return (r=radius, check=check, check_bl=check_bl, coords=coords, sign_r=λ -> -1.0,
        polar=polar,
        domain=(mino=(-Inf, λ_h), endpoint_closed=(false, true),
            endpoint_roles=(:past_repeated_root_asymptote, :future_horizon),
            horizon_lambda=λ_h),
        reference=(lambda0_event=:future_horizon, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=λ_ref, t_phi_zero_radius=rref, tau_zero_event=:future_horizon,
            lambda_regular=λ_h, polar_phase_event=id === :K5 ? :future_horizon : :reference_radius),
        trajectory=(lambda_of_radius=lambda_of_radius,))
end

# ---- assembly -------------------------------------------------------------------------------
function _critical_member(a, energy, lz, q, id, component, classification, polar, phase,
        root, track)
    spec = kerr_geo_case(id)
    equatorial = abs(q) <= 1.0e-12
    name = id === :K1 ? (equatorial ? :isco : :isso) :
        id in CRITICAL_ROLES.on_root ? (equatorial ? :unstable_circular : :unstable_spherical) :
        id === :K4 ? :homoclinic : :whirling
    track = merge(track, (reference=merge(track.reference, (polar_phase=phase,
        polar_phase_convention=polar.metadata.phase_convention)),))
    return _engine_member(:critical, id, a, energy, lz, q, polar, track;
        component=component,
        roots=(repeated=(radius=float(root.Radius), multiplicity=root.Multiplicity,
                derivatives=kerr_radial_derivatives(a, energy, lz, q, root.Radius))),
        status=(name=name, stability=root.Multiplicity == 3 ? :marginal : :unstable,
            formula_family=spec.FormulaFamily, family_member_case_ids=spec.FamilyMemberCaseIds,
            formula_kind=track.model === nothing ? :constant_radius : track.model.kind))
end
