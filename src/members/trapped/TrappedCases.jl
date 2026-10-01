# Class N (Trapped) members N1-N6: classification of E < 0 constants and assembly
# of the full, outgoing and incoming trapped worldlines.

const _CAUSAL_SAMPLES = 257

"""
    KerrGeoTrappedClassification

The result of `kerr_geo_trapped_classify` for E < 0 constants with |a| < 1: the `CaseId`
(N1–N6) and its `DispositionId` (NFD01–NFD06), the `EnergyRegime` and `MetricLimit`, the
radial root structure (`Roots`, as `kerr_geo_root_structure` returns it), the `PolarSector`,
P(r₊) (`HorizonMomentum`), the `TurningRadius`, the Mino
`Domain` between the two horizons, and the `Conditions` and `Status` of the classification.
"""
struct KerrGeoTrappedClassification
    CaseId::Symbol
    DispositionId::Symbol
    EnergyRegime::Symbol
    MetricLimit::Symbol
    Roots::NamedTuple
    PolarSector::Symbol
    HorizonMomentum::Float64
    TurningRadius::Float64
    Domain::NamedTuple
    Conditions::NamedTuple
    Status::NamedTuple
end

# The N case (as its disposition NFD0k) of an E < 0 root structure, or `nothing`.
function _trapped_disposition(energy, structure)
    below_mult = Tuple(item.multiplicity for item in structure.below_horizon)
    exterior_mult = Tuple(item.multiplicity for item in structure.exterior)
    complex_count = structure.complex_root_count

    if -1.0 < energy < 0.0
        below_mult == (1,) && exterior_mult == (1, 1, 1) &&
            complex_count == 0 && return :NFD01
        below_mult == (1,) && exterior_mult == (1, 2) &&
            complex_count == 0 && return :NFD02
        below_mult == (1, 1, 1) && exterior_mult == (1,) &&
            complex_count == 0 && return :NFD03
        below_mult == (1,) && exterior_mult == (1,) &&
            complex_count == 2 && return :NFD04
    elseif energy == -1.0
        below_mult == (1,) && exterior_mult == (1, 1) &&
            complex_count == 0 && return :NFD05
    elseif energy < -1.0
        below_mult == (1, 1) && exterior_mult == (1, 1) &&
            complex_count == 0 && return :NFD06
    end
    return nothing
end

function _polar_zmax2(a, energy, lz, q, sector)
    sector === :equatorial && return 0.0
    roots = _polar_quadratic_roots(a, energy, lz, q)
    iszero(roots.c) && q + lz^2 <= 0 && error("Polar motion with E=-1 is degenerate.")
    roots.disc >= 0 || error("Polar roots with E<0 are not real.")
    return roots.c >= 0 ? roots.u_small : max(roots.u_small, roots.u_big)
end

"""
    kerr_geo_trapped_classify(a, E, Lz, Q)

Classify subextremal E < 0 constants into one of the Class N (Trapped) cases N1–N6, one per
radial root structure; `DispositionId` is the matching label NFD01–NFD06 (Nk = NFD0k).
"""
function kerr_geo_trapped_classify(a::Real, energy::Real, lz::Real, q::Real)
    all(isfinite, (a, energy, lz, q)) || throw(DomainError(
        (a, energy, lz, q), "Trapped constants must be finite."))
    metric = kerr_metric_limit(a)
    metric !== :schwarzschild || throw(DomainError(
        a, "Exterior E<0 timelike motion is absent in Schwarzschild."))
    metric !== :extremal || throw(DomainError(
        a, "Class N at extremal spin |a| = 1 is built by `kerr_geo_extremal`."))
    energy < 0.0 || throw(DomainError(
        energy, "Class N requires E<0."))
    q >= 0.0 || throw(DomainError(
        q, "The admitted Class N polar sector requires Q>=0."))
    a * lz < 0.0 || throw(DomainError(
        lz, "Class N requires angular momentum opposite to the signed spin."))

    horizons = kerr_horizons(a)
    pplus = kerr_radial_momentum(a, energy, lz, horizons.rplus)
    pplus > 0.0 || throw(DomainError(
        pplus, "A future event-horizon crossing requires Pplus>0."))
    polar = kerr_polar_admissibility(a, energy, lz, q)
    polar.admissible || throw(DomainError(
        (energy, lz, q), "The E<0 polar potential is inadmissible."))

    structure = kerr_geo_root_structure(a, energy, lz, q)
    isempty(structure.horizon_coincident) || throw(DomainError(
        structure.horizon_coincident,
        "Class N excludes constants with a radial root on the outer horizon."))
    disposition = _trapped_disposition(float(energy), structure)
    disposition === nothing &&
        error("The constants do not match an admitted Class N radial disposition.")
    case_id = _trapped_case_id(disposition)
    first_exterior = first(structure.exterior)
    first_exterior.multiplicity == 1 || error(
        "The physical Class N turning endpoint must be a simple root.")
    turn = first_exterior.radius

    sector = iszero(q) ? :equatorial : :pendular
    zmax2 = _polar_zmax2(a, energy, lz, q, sector)
    0.0 <= zmax2 < 1.0 + POLAR_ROOT_SLACK || error(
        "The Class N polar turning value lies outside the physical interval.")
    zmax2 = clamp(zmax2, 0.0, 1.0)
    minimum_stationary_limit = 1.0 + sqrt(max(1.0 - a^2 * zmax2, 0.0))
    turn < minimum_stationary_limit || throw(DomainError(
        turn,
        "The full radial component is not strictly confined to the ergoregion for every admitted polar phase."))
    midpoint = 0.5 * (horizons.rplus + turn)
    radial_midpoint = kerr_radial_potential(a, energy, lz, q, midpoint)
    radial_midpoint > 0.0 || error(
        "The horizon-to-turn Class N radial interval is not allowed.")
    outer_probe = turn + max(1.0e-7, 1.0e-5 * max(1.0, abs(turn)))
    radial_outer = kerr_radial_potential(a, energy, lz, q, outer_probe)
    radial_outer < 0.0 || error(
        "The first exterior Class N root does not terminate the allowed interval.")

    domain = (
        radial=(horizons.rplus, turn),
        radial_endpoint_closed=(false, true),
        endpoint_roles=(:past_horizon, :finite_turning_point),
    )
    conditions = (
        energy_below_zero=true,
        nonzero_subextremal_spin=true,
        opposite_signed_angular_momentum=true,
        nonnegative_carter_q=true,
        strict_future_horizon_generator=true,
        polar_admissible=true,
        simple_physical_turn=true,
        radial_interval_allowed=true,
        immediately_outer_interval_forbidden=true,
        all_phase_ergoregion_clearance=minimum_stationary_limit - turn,
        zmax2=zmax2,
    )
    return KerrGeoTrappedClassification(
        case_id,
        disposition,
        kerr_geo_case(case_id).EnergyRegime,
        metric,
        structure,
        sector,
        float(pplus),
        float(turn),
        domain,
        conditions,
        (
            supported=true,
            reason=:trapped_component_classified,
            broad_class=:trapped,
            tier=:primary,
        ),
    )
end


function _build_trapped(a, energy, lz, q, classification;
        component=:full,
        polar_phase::Real=0.0,
        disposition_id=nothing)
    disposition_id === nothing || disposition_id === classification.DispositionId ||
        error("Requested disposition $(disposition_id) does not match $(classification.DispositionId).")
    component in (:full, :outgoing, :incoming) ||
        error("Select component=:full, :outgoing, or :incoming.")
    selected_component = component
    radial = _trapped_radial_model(classification.DispositionId, energy, classification.Roots)
    residues = _radial_residues(a, energy, lz)
    polar = _polar_solution(
        a, energy, lz, q, classification.PolarSector, float(polar_phase))
    lambda_horizon = radial.mino(residues.rplus)
    lambda_horizon > 0.0 || error("Class N requires a positive half-duration.")
    # t, φ, τ (zero at the turning event) and the retarded (u, χ) / advanced (v, ψ) charts
    # (zero on the past / future horizon): radial spectral engine + polar primitive
    rplus = residues.rplus
    radial_radius = radial.radius
    radius_of(lambda) = abs(lambda) >= lambda_horizon ? rplus : radial_radius(abs(lambda))
    radial_potential = _radial_potential_from_roots(a, energy, lz, q, classification.Roots)
    coords = _engine_coordinates(a, energy, lz, q, radius_of, _polar_primitive(polar);
        potential=radial_potential,
        domain=(-lambda_horizon, lambda_horizon), ends=(:horizon, :horizon), turn=0.0,
        σ=-1.0, λ_bl=0.0, λ_regular=lambda_horizon, σ_regular=-1.0)
    retarded = _regular_chart(coords, 1.0, -lambda_horizon)

    function check_full(lambda)
        lam = float(lambda)
        abs(lam) <= lambda_horizon + MINO_ENDPOINT_TOL || throw(DomainError(
            lambda, "Mino time lies outside the full trapped interval."))
        return clamp(lam, -lambda_horizon, lambda_horizon)
    end
    check_full_bl(lambda) = (abs(check_full(lambda)) < lambda_horizon ||
        throw(DomainError(lambda, "BL t and phi exclude both exact horizon endpoints."));
        check_full(lambda))
    function selected(lam, lambda)
        selected_component === :outgoing && lam > MINO_ENDPOINT_TOL && throw(DomainError(
            lambda, "The outgoing component ends at the radial turning event."))
        selected_component === :incoming && lam < -MINO_ENDPOINT_TOL && throw(DomainError(
            lambda, "The incoming component starts at the radial turning event."))
        return selected_component === :outgoing ? min(lam, 0.0) :
            selected_component === :incoming ? max(lam, 0.0) : lam
    end
    check_component(lambda) = selected(check_full(lambda), lambda)
    check_component_bl(lambda) = selected(check_full_bl(lambda), lambda)
    function full_bl(lambda)
        lam = check_full_bl(lambda)
        ps = polar.formula(lam)
        tpt = _coords_tphitau(coords, lam)
        return (
            lambda=lam,
            t=tpt[1],
            r=radius_of(lam),
            theta=ps.theta,
            phi=tpt[2],
            z=ps.z,
            tau=tpt[3],
        )
    end
    function past_regular(lambda)
        lam = check_full(lambda)
        lam <= MINO_ENDPOINT_TOL || throw(DomainError(
            lambda, "The retarded chart (u, χ) covers only the outgoing half, λ ≤ 0."))
        u, chi = retarded(lam)
        return (u=u, chi=chi)
    end
    function future_regular(lambda)
        lam = check_full(lambda)
        lam >= -MINO_ENDPOINT_TOL || throw(DomainError(
            lambda, "The advanced chart (v, ψ) covers only the incoming half, λ ≥ 0."))
        return (v=_coords_v(coords, lam), psi=_coords_psi(coords, lam))
    end
    function outgoing_view(lambda)
        lam = check_full_bl(lambda)
        lam <= MINO_ENDPOINT_TOL || throw(DomainError(
            lambda, "The outgoing view requires lambda<=0."))
        return merge(full_bl(min(lam, 0.0)), past_regular(min(lam, 0.0)))
    end
    function incoming_view(lambda)
        lam = check_full_bl(lambda)
        lam >= -MINO_ENDPOINT_TOL || throw(DomainError(
            lambda, "The incoming view requires lambda>=0."))
        return merge(full_bl(max(lam, 0.0)), future_regular(max(lam, 0.0)))
    end

    # velocities and residuals on the full worldline (λ < 0 outgoing, λ > 0 incoming)
    position = _polar_position(polar)
    kin = _kinematics(a, energy, lz, q, λ -> radius_of(check_full(λ)),
        λ -> radius_of(check_full_bl(λ)), λ -> position(check_full(λ)), λ -> -sign(λ),
        radial_potential)
    k = kin.state
    full_velocity(λ) = (ut=_kin_ut(k, λ), ur=_kin_ur(k, λ), uz=_kin_uz(k, λ), utheta=_kin_utheta(k, λ),
        uphi=_kin_uphi(k, λ), dtau_dlambda=_kin_dtau(k, λ))
    full_residuals(λ) = (radial=_kin_ur(k, λ)^2 - k.R(k.r(λ)),
        polar_z=_kin_uz(k, λ)^2 - kerr_polar_z_potential(a, energy, lz, q, _kin_z(k, λ)),
        normalization=_kin_normalization(k, λ))

    # future-directed and co-rotating everywhere (the ergoregion clearance of the
    # classification guarantees it; checked on a grid of the assembled rates)
    causal_margin = max(1.0e-8, 1.0e-7 * lambda_horizon)
    causal_grid = range(-lambda_horizon + causal_margin, lambda_horizon - causal_margin;
        length=_CAUSAL_SAMPLES)
    minimum_dt = minimum(λ -> _kin_ut(k, λ), causal_grid)
    minimum_corotation = minimum(λ -> a * _kin_uphi(k, λ), causal_grid)
    minimum_dt > 0.0 || error(
        "The Class N worldline is not future-directed (minimum dt/dλ = $(minimum_dt)).")
    minimum_corotation > 0.0 || error(
        "The Class N worldline is not co-rotating (minimum a dφ/dλ = $(minimum_corotation)).")

    full_domain = (-lambda_horizon, lambda_horizon)
    outgoing_domain = (-lambda_horizon, 0.0)
    incoming_domain = (0.0, lambda_horizon)
    selected_domain = selected_component === :full ? full_domain :
        selected_component === :outgoing ? outgoing_domain : incoming_domain
    domain = (
        mino=selected_domain,
        selected_component=selected_component,
        full=full_domain,
        outgoing=outgoing_domain,
        incoming=incoming_domain,
        endpoint_closed=(true, true),
        endpoint_roles=selected_component === :full ? (:past_horizon, :future_horizon) :
            selected_component === :outgoing ? (:past_horizon, :finite_turning_point) :
            (:finite_turning_point, :future_horizon),
        turning_lambda=0.0,
        turning_radius=radial.turn,
        horizon_half_duration=lambda_horizon,
    )
    reference = (
        lambda0_event=:finite_turning_point,
        t_phi_zero_event=:finite_turning_point,
        t_phi_zero_lambda=0.0,
        t_phi_zero_radius=radial.turn,
        tau_zero_event=:finite_turning_point,
        lambda_regular=lambda_horizon,
        lambda_retarded=-lambda_horizon,
        polar_phase=float(polar_phase),
        polar_phase_convention=polar.metadata.phase_convention,
    )
    z(lambda) = position(check_component(lambda))[1]
    trajectory = (
        t=lambda -> _coords_t(coords, check_component_bl(lambda)),
        r=lambda -> radius_of(check_component(lambda)),
        theta=lambda -> acos(clamp(z(lambda), -1.0, 1.0)),
        z=z,
        phi=lambda -> _coords_phi(coords, check_component_bl(lambda)),
        tau=lambda -> _coords_tau(coords, check_component_bl(lambda)),
        u=lambda -> (check_component(lambda); past_regular(lambda).u),
        chi=lambda -> (check_component(lambda); past_regular(lambda).chi),
        v=lambda -> (check_component(lambda); future_regular(lambda).v),
        psi=lambda -> (check_component(lambda); future_regular(lambda).psi),
        lambda_of_radius=radial.mino,
        full=full_bl,
        outgoing=outgoing_view,
        incoming=incoming_view,
    )
    velocity = (
        ut=λ -> _kin_ut(k, check_component_bl(λ)),
        ur=λ -> _kin_ur(k, check_component_bl(λ)),
        uz=λ -> _kin_uz(k, check_component_bl(λ)),
        utheta=λ -> _kin_utheta(k, check_component_bl(λ)),
        uphi=λ -> _kin_uphi(k, check_component_bl(λ)),
        dtau_dlambda=λ -> _kin_dtau(k, check_component_bl(λ)),
        full=full_velocity,
    )
    residuals = (
        radial=λ -> (μ = check_component_bl(λ); _kin_ur(k, μ)^2 - k.R(k.r(μ))),
        polar_z=λ -> (μ = check_component_bl(λ);
            _kin_uz(k, μ)^2 - kerr_polar_z_potential(a, energy, lz, q, _kin_z(k, μ))),
        normalization=λ -> _kin_normalization(k, check_component_bl(λ)),
        full=full_residuals,
    )
    status = (
        supported=true,
        disposition_id=classification.DispositionId,
        energy_regime=classification.EnergyRegime,
        energy_sign=kerr_energy_sign(energy),
        metric_limit=classification.MetricLimit,
        selected_component=selected_component,
        formula_kind=radial.kind,
        # past half: retarded (u, χ), zero on the white-hole horizon; future half: advanced
        # (v, ψ), zero on the black-hole horizon; t → ∓∞ and φ diverges at both ends
        charts=(past=:retarded_u_chi, future=:advanced_v_psi),
        causality=(
            samples=_CAUSAL_SAMPLES,
            minimum_dt_dlambda=minimum_dt,
            minimum_a_dphi_dlambda=minimum_corotation,
            all_phase_ergoregion_clearance=
                classification.Conditions.all_phase_ergoregion_clearance,
        ),
        polar=polar.metadata,
        classification=classification,
    )
    return _member(:trapped, classification.CaseId, kerr_geo_tier(classification.CaseId),
        nothing, (a=float(a), E=float(energy), Lz=float(lz), Q=float(q)),
        (radial=radial.roots, polar=polar.metadata), reference, domain, trajectory, velocity,
        (radial=radial_potential, polar_z=kin.Θ), residuals, status,
        SpectralStatus(() -> (_coords_spectral(coords), _polar_spectral(polar))))
end

"""
    kerr_geo_trapped(a, E, Lz, Q; component=:full, polar_phase=0.0, disposition_id=nothing)

The Trapped member (cases N1–N6) of E < 0 constants with |a| < 1. The worldline leaves the
past (white-hole) horizon at `λ = −λH`, turns at λ = 0 (where t, φ and τ vanish) and crosses
the future (black-hole) horizon at `λ = λH`, inside the ergoregion throughout. `component` selects the full
worldline (`:full`), the outgoing half (`:outgoing`, λ ≤ 0) or the incoming half (`:incoming`,
λ ≥ 0); (u, χ) vanish on the past horizon and (v, ψ) on the future horizon. A given
`disposition_id` (NFD01–NFD06) must match the classification. `Status.causality` records the
minima of dt/dλ and a dφ/dλ over the worldline (both positive: future-directed and co-rotating).
"""
function kerr_geo_trapped(a::Real, energy::Real, lz::Real, q::Real; kwargs...)
    classification = kerr_geo_trapped_classify(a, energy, lz, q)
    return _build_trapped(a, energy, lz, q, classification; kwargs...)
end

function kerr_geo_trapped(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_trapped(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_trapped(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_trapped(a, constants...; kwargs...)
end

"""
    kerr_geo_trapped_case(case_id, a, E, Lz, Q; kwargs...)

`kerr_geo_trapped(a, E, Lz, Q; kwargs...)`, checked to be the Trapped case `case_id`
(`:N1`–`:N6`).
"""
function kerr_geo_trapped_case(
        case_id::Symbol, a::Real, energy::Real, lz::Real, q::Real; kwargs...)
    case_id in TRAPPED_CASE_IDS || error("Class N cases are N1-N6, not $(case_id).")
    orbit = kerr_geo_trapped(a, energy, lz, q; kwargs...)
    orbit.CaseId === case_id || error("These constants are case $(orbit.CaseId), not $(case_id).")
    return orbit
end

function kerr_geo_trapped_case(
        case_id::Symbol, a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_trapped_case(
        case_id, a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_trapped_case(case_id::Symbol, a::Real,
        constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_trapped_case(case_id, a, constants...; kwargs...)
end
