# Family assembly: `kerr_geodesic` classifies the constants and builds every member they
# admit, one slot per broad class (Stable, Critical, Plunge, Capture, Scatter, Trapped).

"""
    KerrGeodesicFamily

Output of `kerr_geodesic`: every member admitted by one set of constants, in the slot of its
broad class. `Stable`, `Plunge`, `Capture`, `Scatter` and `Trapped` hold one member or
`nothing`; `Critical` is a tuple of the Critical members ordered by role (on the root, outer
side, inner side). Horizon-root (H) and extremal (X) members sit in the slot of their
class. `BroadClass` summarises the family (`:critical` when it has Critical members).
Construction failures are recorded in `Status.member_errors` when present.
The affected slot is empty; a returned family need not contain every admitted member.
"""
struct KerrGeodesicFamily
    InputType::Symbol
    Parameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    RootClass::String
    Stable::Union{Nothing,KerrGeoStableComponent}
    Critical::Tuple{Vararg{KerrGeoCriticalComponent}}
    Plunge::Union{Nothing,KerrGeoPlungeComponent}
    Capture::Union{Nothing,KerrGeoCaptureComponent}
    Scatter::Union{Nothing,KerrGeoScatterComponent}
    Trapped::Union{Nothing,KerrGeoTrappedComponent}
    BroadClass::Symbol
    Status::NamedTuple
end

"""
    kerr_geo_members(kg)

The members of the family `kg` as a tuple, in class order: Stable, the Critical members,
Plunge, Capture, Scatter, Trapped.
"""
function kerr_geo_members(f::KerrGeodesicFamily)
    members = Any[]
    for row in KERR_GEO_CLASSES
        slot = getfield(f, row.slot)
        slot isa Tuple ? append!(members, slot) : slot === nothing || push!(members, slot)
    end
    return Tuple(members)
end


_family(input, parameters, constants, root_class, broad, status; stable=nothing,
        critical=(), plunge=nothing, capture=nothing, scatter=nothing, trapped=nothing) =
    KerrGeodesicFamily(input, parameters, constants, root_class, stable, Tuple(critical),
        plunge, capture, scatter, trapped, broad, status)

function Base.show(io::IO, kg::KerrGeodesicFamily)
    print(io, "KerrGeodesicFamily(", kg.BroadClass, ", cases=")
    show(io, get(kg.Status, :case_ids, ()))
    print(io, ", member_errors=", length(get(kg.Status, :member_errors, ())), ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeodesicFamily)
    println(io, "KerrGeodesicFamily (", kg.InputType, ")")
    _show_summary_field(io, "Constants", kg.ConstantsOfMotion)
    _show_summary_field(io, "Parameters", kg.Parameters)
    for row in KERR_GEO_CLASSES
        slot = getfield(kg, row.slot)
        if slot isa Tuple
            isempty(slot) || _show_summary_field(io, row.name, Tuple(m.CaseId for m in slot))
        elseif slot !== nothing
            _show_summary_field(io, row.name, slot.CaseId)
        end
    end
    _show_summary_status(io, kg.Status)
end

"""
    kerr_geodesic(a, (E, Lz, Q); kwargs...)
    kerr_geodesic(a; constants=(E, Lz, Q), kwargs...)
    kerr_geodesic(a, p, e, x; input=:apex, kwargs...)

Build the `KerrGeodesicFamily` of spin `a` and constants of motion `(E, Lz, Q)` (a tuple, or
a NamedTuple with fields `E`, `Lz`, `Q`): every member these constants admit, each in the
slot of its broad class. The constants are used exactly as given.

The four-argument form takes APEX parameters (semi-latus rectum ``p``, eccentricity ``e``,
``x = \\cos\\iota``) and converts them to `(E, Lz, Q)`; a spin within `8eps()` of ``\\pm 1`` is taken as
exactly ``\\pm 1``, and the original input is kept in `Status.input_provenance`. A bound eccentric
orbit keeps the turning points ``p/(1 \\mp e)`` in its Stable member; `Status.apex_root_geometry`
records whether they were used (`accepted`) or why not (`reason`), and
`Status.component_root_models` the roots each component is built from. With
`input=:constants` the three numbers are read as `(E, Lz, Q)`.

Returns a `KerrGeodesicFamily`. Each member (a `KerrGeoComponent`) gives the Boyer–Lindquist
coordinates as functions of Mino time ``\\lambda``, defined by ``d\\tau/d\\lambda = \\Sigma = r^2 + a^2\\cos^2\\theta``:
`m.Trajectory.t(λ)`, `r`, `theta`, `phi`, `tau`, and the rates ``dx^\\mu/d\\lambda`` in `m.Velocity`.

```julia
kg = kerr_geodesic(0.9, (0.94, 0.1, 12.0))
kg.Plunge.Trajectory.r(0.5)
```

Keywords: `case_id`, `initial_radius`, `radial_sign`, `endpoint_intent` select one member
(exactly one must match); `polar_sector`, `polar_phase`, `polar_hemisphere` fix the polar
motion; `reference_radius` places the zero of ``t`` and ``\\phi`` on Critical and Capture members;
`axis` (`:north`, `:south`) selects motion along the spin axis, at constant azimuth `phi0`;
`initPhases` sets the phases of the Stable member; `trapped_component` (`:full`,
`:outgoing`, `:incoming`) selects the part of the Trapped member, and `disposition_id` is
checked against it.

The orbits are computed in the floating-point type of the input (`Float64`, or `BigFloat` at
the precision of the given numbers); `precision = p` converts the input to `BigFloat` of `p`
bits. The members' functions evaluate at the precision they were built with.
"""
function kerr_geodesic(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geodesic(a, (constants.E, constants.Lz, constants.Q); kwargs...)
end

# The Plunge, Capture or Scatter member of classified constants, or `nothing`.
function _family_member(class, a, energy, lz, q, classification, kwargs)
    ids = [c.CaseId for c in classification.Components
        if c.BroadClass === class && c.CaseId !== nothing]
    isempty(ids) && return nothing
    length(ids) == 1 || error("Multiple $(kerr_geo_class(class).name) components require " *
        "an explicit selection (case_id).")
    id = only(ids)
    sector = get(kwargs, :polar_sector, nothing)
    phase = get(kwargs, :polar_phase, 0.0)
    hemisphere = get(kwargs, :polar_hemisphere, :north)
    reference_radius = get(kwargs, :reference_radius, nothing)
    class === :plunge && return _plunge_component(a, energy, lz, q, classification, id;
        polar_phase=phase, reference_radius=reference_radius)
    class === :scatter && return _scatter_component(a, energy, lz, q, classification, id;
        polar_sector=sector, polar_phase=phase, polar_hemisphere=hemisphere)
    return _capture_component(a, energy, lz, q, classification, id; polar_sector=sector,
        polar_phase=phase, polar_hemisphere=hemisphere, reference_radius=reference_radius)
end

# Motion along the spin axis (Lz = 0, Q = a²(1 - E²)) with an `axis` initial condition: the
# Plunge (E < 1) or Capture (E ≥ 1) axis-infall member. (The classification was made with
# `polar_sector = :axis_constant`, which holds only for axis constants.)
function _family_axis_member(a, energy, classification, kwargs)
    classification.PolarMetadata.selected === :axis_constant || error(
        "Axis initial data requires the axis_constant polar sector.")
    constants = classification.ConstantsOfMotion
    iszero(constants.Lz) && constants.Q == kerr_axis_carter_q(a, energy) || error(
        "Axis construction must preserve the supplied Lz and Q.")
    component = _axis_component(classification, energy < 1 ? :plunge : :capture)
    return _axis_infall_member(a, energy, component; axis=kwargs[:axis],
        phi0=get(kwargs, :phi0, 0.0), reference_radius=get(kwargs, :reference_radius, nothing))
end

# The Critical members of the classified components, in role order.
function _family_critical_members(a, energy, lz, q, classification, kwargs, errors)
    members = KerrGeoCriticalComponent[]
    for component in classification.Components
        # an unassigned Critical interval (no case ID) stays visible in Status.components
        component.BroadClass === :critical && component.CaseId !== nothing || continue
        member = _isolated_member(() -> _critical_component(a, energy, lz, q,
                classification, component.CaseId;
                polar_sector=get(kwargs, :polar_sector, nothing),
                polar_phase=get(kwargs, :polar_phase, 0.0),
                polar_hemisphere=get(kwargs, :polar_hemisphere, :north),
                reference_radius=get(kwargs, :reference_radius, nothing)),
            component.CaseId, errors, kwargs)
        member === nothing || push!(members, member)
    end
    return Tuple(members)
end

# One-word summary of a family, from its case ids (falls back to the outcome at infinity);
# `:none` when the constants admit no member.
function _family_broad_class(ids, outcome)
    isempty(ids) && return :none
    any(id -> kerr_geo_case_class(id) === :trapped, ids) && return :trapped
    any(id -> kerr_geo_case_class(id) === :critical, ids) && return :critical
    any(id -> id in (:D1, :D2), ids) && return :scatter
    any(id -> id in (:B3, :B4), ids) && return :plunge
    any(id -> id in (:C1, :C2, :C3, :C4), ids) && return :capture
    return outcome.outcome
end

function _family_selection(classification, kwargs)
    case_id = get(kwargs, :case_id, nothing)
    initial_radius = get(kwargs, :initial_radius, nothing)
    radial_sign = get(kwargs, :radial_sign, nothing)
    endpoint_intent = get(kwargs, :endpoint_intent, nothing)
    requested = any(value -> value !== nothing,
        (case_id, initial_radius, radial_sign, endpoint_intent))
    requested || return (
        requested=false,
        component=nothing,
        case_id=nothing,
        criteria=(case_id=nothing, initial_radius=nothing,
            radial_sign=nothing, endpoint_intent=nothing),
    )
    component = kerr_geo_select_component(
        classification.Components;
        case_id=case_id,
        initial_radius=initial_radius,
        radial_sign=radial_sign,
        endpoint_intent=endpoint_intent,
    )
    return (
        requested=true,
        component=component,
        case_id=component.CaseId,
        criteria=(case_id=case_id, initial_radius=initial_radius,
            radial_sign=radial_sign, endpoint_intent=endpoint_intent),
    )
end

function _exact_extremal_family(a, energy, lz, q, kwargs)
    polar_sector = get(kwargs, :polar_sector, nothing)
    polar_phase = get(kwargs, :polar_phase, 0.0)
    axis = get(kwargs, :axis, nothing)
    if axis !== nothing
        polar_sector in (nothing, :axis_constant) || error(
            "Axis initial data conflicts with the requested polar sector $(polar_sector).")
        iszero(lz) && q == kerr_axis_carter_q(a, energy) || error(
            "Axis construction must preserve the supplied Lz and Q.")
    end
    reference_radius = get(kwargs, :reference_radius, nothing)
    exact = kerr_geo_extremal_family(
        a, energy, lz, q;
        polar_sector=polar_sector,
        polar_phase=polar_phase,
        polar_hemisphere=get(kwargs, :polar_hemisphere, :north),
        axis=axis,
        reference_radius=reference_radius,
    )
    members = exact.Members
    of_class(c) = Tuple(member for member in members if kerr_geo_member_class(member) === c)
    single(c) = (found = of_class(c); length(found) <= 1 || error(
        "The exact-extremal family contains several $(c) members."); isempty(found) ? nothing : only(found))
    case_ids = exact.Classification.case_ids
    broad_class = _family_broad_class(case_ids, _outcome_at_infinity(a, energy, lz, q))
    metric_limit = a == 1 ? :extremal_plus : :extremal_minus
    requested_case = get(kwargs, :case_id, nothing)
    selected_member = if requested_case === nothing
        length(members) == 1 ? only(members) : nothing
    else
        matches = filter(member -> member.CaseId === requested_case, members)
        length(matches) == 1 || error(
            "Case $(requested_case) is not an admitted exact-extremal member.")
        only(matches)
    end
    status = (
        supported=!isempty(members),
        reason=isempty(members) ?
            (a == 1 ? :no_admitted_exact_positive_component :
             :no_admitted_exact_negative_component) :
            (a == 1 ? :exact_positive_components_available :
             :exact_negative_components_available),
        classification=exact.Classification,
        case_ids=case_ids,
        components=members,
        selected_case=selected_member === nothing ? nothing : selected_member.CaseId,
        selected_component=selected_member,
        selection_hint=selected_member !== nothing ?
            (case_id=selected_member.CaseId, stage=:optional_selection) :
            length(case_ids) <= 1 ? nothing :
            (case_ids=case_ids,
             instruction=:select_from_family_members_by_case_id),
        metric_limit=metric_limit,
    )
    return _family(:constants, (a=a,), (E=energy, Lz=lz, Q=q),
        isempty(case_ids) ? "ExactExtremalExcluded" :
            (a == 1 ? "ExactExtremalPlus[" : "ExactExtremalMinus[") *
                join(String.(case_ids), ",") * "]",
        broad_class, status;
        stable=single(:stable), critical=of_class(:critical), plunge=single(:plunge),
        capture=single(:capture), scatter=single(:scatter), trapped=single(:trapped))
end

# The constants are used exactly as given: the quartic coefficient E² − 1 decides whether
# an orbit reaches infinity however small it is, so no energy is moved to E = 1. Only the
# APEX conversion rounds a spin within a few ulps of |a| = 1 to exactly ±1 (below).
function kerr_geodesic(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; precision=nothing,
        kwargs...)
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geodesic(BigFloat(a), BigFloat.(constants); kwargs...)
    end
    T = _float_type(a, constants...)
    return _with_precision(T, _input_precision(a, constants...)) do
        _kerr_geodesic(T(a), T.(constants); kwargs...)
    end
end

# `stable_geometry`: the APEX turning-point root geometry the Stable member is built from
# (`_apex_root_geometry`); every other member uses the roots of the constants.
function _kerr_geodesic(a, constants; stable_geometry=nothing, kwargs...)
    energy, lz, q = constants
    if haskey(kwargs, :axis)
        requested = get(kwargs, :polar_sector, nothing)
        sector = _polar_sector(kerr_polar_sector_candidates(a, energy, lz, q),
            requested === nothing ? :axis_constant : requested)
        sector === :axis_constant || error(
            "Axis initial data conflicts with the requested polar sector $(sector).")
        iszero(lz) && q == kerr_axis_carter_q(a, energy) || error(
            "Axis construction must preserve the supplied Lz and Q.")
    end
    family = abs(a) == 1 ? _exact_extremal_family(a, energy, lz, q, kwargs) :
        energy < 0 ? _trapped_family(a, energy, lz, q, kwargs) :
        _horizon_root_family(a, energy, lz, q, kwargs)
    family === nothing && (family = _classified_family(a, energy, lz, q, kwargs;
        stable_geometry=stable_geometry))
    _validate_selection(family, kwargs)
    return family
end

# Selection keywords (case_id, initial_radius, radial_sign, endpoint_intent) are checked the
# same way on every tier: exactly one built member must match them, and it must be the
# family's selected member.
function _validate_selection(family, kwargs)
    criteria = (case_id=get(kwargs, :case_id, nothing),
        initial_radius=get(kwargs, :initial_radius, nothing),
        radial_sign=get(kwargs, :radial_sign, nothing),
        endpoint_intent=get(kwargs, :endpoint_intent, nothing))
    all(isnothing, values(criteria)) && return nothing
    members = collect(kerr_geo_members(family))
    views = [_selection_view(m) for m in members]
    isempty(views) && error("No member of these constants matches the selection $(criteria).")
    chosen = members[findfirst(v -> v === kerr_geo_select_component(views; criteria...), views)]
    selected = get(family.Status, :selected_case, nothing)
    selected === nothing || selected === chosen.CaseId || error(
        "The selection $(criteria) matches $(chosen.CaseId), not the family's $(selected).")
    return chosen
end

# The selection record of a member: its classified component, or (horizon- and extremal-tier
# members, which carry none) the same fields derived from its domain and roots.
function _selection_view(m)
    m.Component === nothing || return m.Component
    lo, hi = float.(m.Domain.mino)
    r = m.Trajectory.r
    λmid = isfinite(lo) && isfinite(hi) ? (lo + hi) / 2 : isfinite(lo) ? lo + 1 :
        isfinite(hi) ? hi - 1 : 0.0
    rmid = r(λmid)
    roots = sort!([x isa Number ? float(x) : float(x.radius) for x in get(m.Roots, :radial, ())])
    rplus = _rplus(m.ConstantsOfMotion.a)
    at_infinity(role) = role in (:past_infinity, :future_infinity)
    at_horizon(role) = role in (:past_horizon, :future_horizon,
        :past_horizon_root_asymptote, :future_horizon_root_asymptote)
    at_repeated_root(role) = role in (:past_repeated_root_asymptote,
        :future_repeated_root_asymptote)
    at_asymptote(role) = at_repeated_root(role) || role in
        (:past_horizon_root_asymptote, :future_horizon_root_asymptote)
    function endpoint(role, λ, side)
        at_infinity(role) && return (Radius=Inf, Included=false, Kind=:infinity)
        role in (:past_horizon, :future_horizon) &&
            return (Radius=rplus, Included=true, Kind=:outer_horizon)
        isfinite(λ) && return (Radius=r(λ), Included=true, Kind=:radial_root)
        near = side < 0 ? filter(x -> x <= rmid, roots) : filter(x -> x >= rmid, roots)
        return (Radius=isempty(near) ? rmid : (side < 0 ? maximum(near) : minimum(near)),
            Included=false, Kind=:radial_root)
    end
    past, future = m.Domain.endpoint_roles
    e1 = endpoint(past, lo, -1); e2 = endpoint(future, hi, +1)
    lower, upper = e1.Radius <= e2.Radius ? (e1, e2) : (e2, e1)
    orientation = m.Role === :on_root || lower.Radius == upper.Radius ? :constant_radius :
        past in (:infinite_past_worldline, :infinite_future_worldline) ? :libration :
        at_infinity(past) && at_infinity(future) ? :inbound_turn_outbound :
        at_horizon(past) && at_horizon(future) ? :inbound_turn_outbound :
        at_repeated_root(past) && at_repeated_root(future) ? :libration :
        at_horizon(past) ? :outward :
        at_asymptote(future) || at_asymptote(past) ? :inward_asymptotic : :inward
    return (CaseId=m.CaseId, BroadClass=kerr_geo_member_class(m), LowerEndpoint=lower,
        UpperEndpoint=upper, RadialOrientation=orientation)
end

# E < 0: the only admitted cases are the Trapped N1-N6 (one per radial root structure).
# Constants that admit no future-directed exterior motion give an empty family.
function _trapped_family(a, energy, lz, q, kwargs)
    requested_case = get(kwargs, :case_id, nothing)
    classification = try
        kerr_geo_trapped_classify(a, energy, lz, q)
    catch err
        err isa DomainError || rethrow()
        return _empty_family(a, energy, lz, q, err.msg)
    end
    trapped = _build_trapped(a, energy, lz, q, classification;
        component=get(kwargs, :trapped_component, :full),
        polar_phase=get(kwargs, :polar_phase, 0.0),
        disposition_id=get(kwargs, :disposition_id, nothing))
    case_id = trapped.CaseId
    requested_case === nothing || requested_case === case_id || error(
        "E < 0 constants admit $(case_id), not $(requested_case)" *
        ".")
    requested_polar = get(kwargs, :polar_sector, nothing)
    requested_polar === nothing ||
        requested_polar === trapped.Status.classification.PolarSector || error(
            "Requested polar sector does not match the $(case_id) constants.")
    status = (
        supported=true,
        reason=:trapped_component_available,
        classification=trapped.Status.classification,
        case_ids=(case_id,),
        components=(trapped,),
        selected_case=case_id,
        selected_component=trapped,
        selection_hint=(
            case_id=case_id,
            component=trapped.Status.selected_component,
            polar_sector=trapped.Status.classification.PolarSector,
        ),
    )
    return _family(:constants, (a=a,), (E=energy, Lz=lz, Q=q),
        String(trapped.Status.disposition_id), :trapped, status; trapped=trapped)
end

# A turning point on the outer horizon (P(r+) = 0, sub-extremal spin, E > 0) has its own
# horizon-root stable or scatter member. Returns `nothing` otherwise.
function _horizon_root_family(a, energy, lz, q, kwargs)
    metric_limit = kerr_metric_limit(a)
    if 0 < abs(a) < 1 && energy > 0 &&
            _horizon_root(a, energy, lz)
        horizons = kerr_horizons(a)
        pplus = kerr_radial_momentum(a, energy, lz, horizons.rplus)
        horizon_structure = kerr_geo_root_structure(a, energy, lz, q)
        if !isempty(horizon_structure.horizon_coincident)
            polar_sector = get(kwargs, :polar_sector, nothing)
            polar_phase = get(kwargs, :polar_phase, 0.0)
            member = try
                if energy < 1
                    _horizon_stable(a, energy, lz, q, horizons, horizon_structure;
                        polar_sector=polar_sector, polar_phase=polar_phase)
                else
                    _horizon_scatter(a, energy, lz, q, horizon_structure;
                        polar_sector=polar_sector, polar_phase=polar_phase,
                        polar_hemisphere=get(kwargs, :polar_hemisphere, :north))
                end
            catch err
                err isa Union{ErrorException,DomainError} || rethrow()
                return _empty_family(a, energy, lz, q,
                    _horizon_coincident_reason(a, energy, lz, q, horizons.rplus,
                        sprint(showerror, err)))
            end
            is_stable = member isa KerrGeoStableComponent
            is_scatter = member isa KerrGeoScatterComponent
            status = (
                supported=true,
                reason=is_stable ? :horizon_root_stable_component_available :
                    :horizon_root_scatter_component_available,
                case_ids=(member.CaseId,),
                components=(member,),
                selected_case=member.CaseId,
                selected_component=member,
                metric_limit=metric_limit,
                horizon_momentum=pplus,
                root_structure=horizon_structure,
            )
            return _family(:constants, (a=a,), (E=energy, Lz=lz, Q=q),
                String(member.CaseId), is_stable ? :stable : :scatter, status;
                stable=is_stable ? member : nothing,
                scatter=is_scatter ? member : nothing)
        end
    end
    return nothing
end

# P(r+) = 0 constants that admit no horizon-tier member. If R > 0 just outside r+
# the allowed region reaches the horizon at a turning point: that orbit runs through the
# bifurcation sphere (finite Boyer–Lindquist t at r = r+), which is not modelled here.
function _horizon_coincident_reason(a, energy, lz, q, rplus, message)
    touches = kerr_radial_potential(a, energy, lz, q, rplus * (1 + 1e-6)) > 0
    return touches ?
        "P(r+) = 0 and the allowed region reaches the horizon at a turning point: the " *
        "orbit passes through the bifurcation sphere, which is not modelled ($message)" :
        "P(r+) = 0 and no exterior component is admitted ($message)"
end

_empty_family(a, energy, lz, q, reason) = _family(:constants, (a=a,),
    (E=energy, Lz=lz, Q=q), "none", :none,
    (supported=false, reason=reason, case_ids=(), components=()))

# Generic constants: classify the radial components (A/B/C/D cases) and build every
# member of the family.
# A member that cannot be built does not take the other members of the family with it: the
# failure is recorded in Status.member_errors and that slot stays empty. Errors caused by an
# explicit user selection (case_id, polar_sector, ...) are thrown.
function _isolated_member(build, label, errors, kwargs)
    explicit = any(haskey(kwargs, k) for k in (:case_id, :polar_sector,
        :initial_radius, :radial_sign))
    explicit && return build()
    try
        return build()
    catch err
        (err isa ErrorException || err isa DomainError || err isa ArgumentError) || rethrow()
        push!(errors, (member=label, message=first(split(sprint(showerror, err), '\n'))))
        return nothing
    end
end

function _classified_family(a, energy, lz, q, kwargs; stable_geometry=nothing)
    member_errors = NamedTuple[]
    # Motion exactly along the spin axis (Lz = 0, Q = a²(1 - E²); for a = 0 purely radial
    # motion) defaults to the northern axis. A defaulted axis member that cannot be built is
    # recorded like any other member; a requested one (`axis` given) raises.
    axis_requested = haskey(kwargs, :axis)
    if !axis_requested &&
            kerr_polar_sector_candidates(a, energy, lz, q) == (:axis_constant,)
        kwargs = pairs(merge(NamedTuple(kwargs), (axis=:north,)))
    end
    # a requested sector is checked against the constants even where the axis is the default
    requested_polar = get(kwargs, :polar_sector, haskey(kwargs, :axis) ? :axis_constant : nothing)
    classification = kerr_geo_classify(
        a, energy, lz, q; polar_sector=requested_polar)
    # the Stable member's own classification: of the APEX turning-point geometry if given
    stable_classification = stable_geometry === nothing ? classification :
        kerr_geo_classify(a, energy, lz, q; polar_sector=requested_polar,
            structure=stable_geometry.structure)
    stable_ids = (:A1, :A2)
    components = stable_geometry === nothing ? classification.Components :
        [filter(c -> c.CaseId in stable_ids, stable_classification.Components);
         filter(c -> !(c.CaseId in stable_ids), classification.Components)]
    component_models = Tuple((case_id=c.CaseId, roots=stable_geometry !== nothing &&
        c.CaseId in stable_ids ? :apex_turning_points : :constants) for c in components)
    axis_member = !haskey(kwargs, :axis) ? nothing : axis_requested ?
        _family_axis_member(a, energy, classification, kwargs) :
        _isolated_member(() -> _family_axis_member(a, energy, classification, kwargs),
            energy < 1 ? :plunge : :capture, member_errors, kwargs)
    critical = _family_critical_members(a, energy, lz, q, classification, kwargs, member_errors)
    plunge = axis_member !== nothing && energy < 1 ? axis_member :
        _isolated_member(() -> _family_member(:plunge, a, energy, lz, q, classification, kwargs),
            :plunge, member_errors, kwargs)
    capture = axis_member !== nothing && energy >= 1 ? axis_member :
        _isolated_member(() -> _family_member(:capture, a, energy, lz, q, classification, kwargs),
            :capture, member_errors, kwargs)
    scatter = axis_member === nothing ?
        _isolated_member(() -> _family_member(:scatter, a, energy, lz, q, classification, kwargs),
            :scatter, member_errors, kwargs) : nothing
    stable = any(component -> component.CaseId in stable_ids, stable_classification.Components) ?
        _isolated_member(() -> _stable_component(a, energy, lz, q, stable_classification, nothing;
            initPhases=get(kwargs, :initPhases, (0.0, 0.0, 0.0, 0.0))), :stable,
            member_errors, kwargs) : nothing
    outcome = _outcome_at_infinity(a, energy, lz, q)
    selection = _family_selection((Components=components,), kwargs)
    members = filter(!isnothing, [stable, critical..., plunge, capture, scatter])
    status = (
        supported=any(m -> m.Status.supported, members),
        reason=isempty(members) ? outcome.reason : :components_available,
        classification=classification,
        case_ids=Tuple(c.CaseId for c in components if c.CaseId !== nothing),
        components=components,
        component_root_models=component_models,
        selected_case=selection.case_id,
        selected_component=selection.component,
        selection_hint=selection.requested ?
            (criteria=selection.criteria, case_id=selection.case_id,
             stage=:optional_selection) : classification.SelectionHint,
        member_errors=Tuple(member_errors),
    )
    return _family(:constants, (a=a,), (E=energy, Lz=lz, Q=q), String(outcome.root_class),
        _family_broad_class(classification.CaseIds, outcome), status;
        stable=stable, critical=critical, plunge=plunge, capture=capture, scatter=scatter)
end

function kerr_geodesic(a::Real, p::Real, e::Real, x::Real; input::Symbol=:apex,
        precision=nothing, kwargs...)
    if input == :constants
        return kerr_geodesic(a, (p, e, x); precision=precision, kwargs...)
    elseif input != :apex
        error("Unknown input type. Use :apex or :constants.")
    end
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geodesic(BigFloat(a), BigFloat(p), BigFloat(e), BigFloat(x); kwargs...)
    end
    T = _float_type(a, p, e, x)
    return _with_precision(T, _input_precision(a, p, e, x)) do
        _kerr_geodesic_apex(T(a), T(p), T(e), T(x); kwargs...)
    end
end

function _kerr_geodesic_apex(a, p, e, x; kwargs...)

    # APEX input only: a spin a few ulps from ±1 (a conversion artefact) is the extremal one
    a_input = a
    abs(abs(a) - 1) <= _spin_snap_tol(_float_type(a)) && (a = copysign(one(float(a)), a))
    constants = kerr_geo_constants_of_motion(a, p, e, x)
    energy = constants["E"]
    lz = constants["Lz"]
    q = constants["Q"]
    constants_tuple = (E=energy, Lz=lz, Q=q)
    # The turning points are part of the input: the Stable member keeps them (with the inner
    # roots of the same constants) unless the geometry is not resolved at this precision.
    geometry = _apex_root_geometry(a, p, e, x, energy, lz, q)
    family = _kerr_geodesic(a, (energy, lz, q);
        stable_geometry=geometry.accepted ? geometry : nothing, kwargs...)
    status = merge(family.Status, (
        apex_root_geometry=(accepted=geometry.accepted, reason=geometry.reason,
            diagnostics=geometry.diagnostics),
        input_provenance=(
            kind=:apex,
            original=(a=a_input, p=p, e=e, x=x),
            adapter=:apex_to_canonical_constants,
        ),
        selection_hint=family.Status.selection_hint,
    ))
    return _family(:apex, (a=a, p=p, e=e, x=x), constants_tuple, family.RootClass,
        family.BroadClass, status; stable=family.Stable, critical=family.Critical,
        plunge=family.Plunge, capture=family.Capture, scatter=family.Scatter,
        trapped=family.Trapped)
end

function kerr_geodesic(a::Real; constants=nothing, kwargs...)
    constants === nothing && error("Provide constants=(E,Lz,Q) for one-argument kerr_geodesic.")
    return kerr_geodesic(a, constants; kwargs...)
end
