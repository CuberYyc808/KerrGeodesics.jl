using Base64
using Printf

const KG56_NOTEBOOK_DIR = @__DIR__
const KG56_PACKAGE_ROOT = normpath(joinpath(KG56_NOTEBOOK_DIR, ".."))
# root structure, allowed interval and formula family of each catalogue case
const KG56_REGISTRY_PATH = joinpath(KG56_NOTEBOOK_DIR, "data", "catalogue_registry.tsv")

KG56_PACKAGE_ROOT in LOAD_PATH || pushfirst!(LOAD_PATH, KG56_PACKAGE_ROOT)
if !isdefined(Main, :KerrGeodesics)
    include(joinpath(KG56_PACKAGE_ROOT, "src", "KerrGeodesics.jl"))
end
using .KerrGeodesics

function _kg56_read_registry(path=KG56_REGISTRY_PATH)
    lines = filter(line -> !startswith(line, "#"), readlines(path))
    header = Symbol.(split(first(lines), '\t'; keepempty=true))
    rows = Dict{Symbol,NamedTuple}()
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        values = split(line, '\t'; keepempty=true)
        row = NamedTuple{Tuple(header)}(Tuple(values))
        rows[kerr_geo_case_symbol(row.case)] = row
    end
    return rows
end

const KG56_REGISTRY = _kg56_read_registry()

# One short motion label per catalogue case (plot titles and gallery captions).
const KG56_MOTION_LABELS = Dict{Symbol,String}(
    :A1 => "stable libration",
    :A2 => "stable spherical",
    :K1 => "ISSO",
    :K3 => "unstable spherical",
    :K4 => "homoclinic orbit",
    :K6 => "parabolic unstable spherical",
    :K9 => "hyperbolic unstable spherical",
    :A_H1 => "detached libration",
    :A_H2 => "detached stable spherical",
    :A_X1 => "extremal stable island",
    :A_X2 => "extremal stable spherical",
    :B1 => "inner-root plunge",
    :B2 => "stable-root plunge",
    :K2 => "from the triple root into the horizon",
    :K5 => "whirl from the unstable root into the horizon",
    :B3 => "outer-root plunge",
    :B7 => "lower-double plunge",
    :B8 => "middle-double plunge",
    :B9 => "triple-inner plunge",
    :B4 => "complex-pair plunge",
    :K8 => "parabolic whirl into the horizon",
    :B5 => "parabolic plunge",
    :K11 => "hyperbolic whirl into the horizon",
    :B6 => "hyperbolic plunge",
    :N1 => "trapped plunge",
    :N2 => "trapped double-root plunge",
    :N3 => "trapped outer-turn plunge",
    :N4 => "trapped complex-root plunge",
    :N5 => "trapped cubic plunge",
    :N6 => "trapped hyperbolic plunge",
    :B_X1 => "double-horizon plunge",
    :B_X2 => "triple-horizon plunge",
    :C1 => "parabolic complex-root capture",
    :C2 => "parabolic real-root capture",
    :C6 => "parabolic lower-double capture",
    :C7 => "parabolic upper-double capture",
    :C8 => "parabolic triple-root capture",
    :K7 => "parabolic whirl in from infinity",
    :C5 => "hyperbolic four-complex capture",
    :C3 => "hyperbolic mixed-root capture",
    :C11 => "hyperbolic double-root capture",
    :C4 => "hyperbolic four-real capture",
    :C9 => "hyperbolic middle-double capture",
    :C10 => "hyperbolic upper-double capture",
    :C12 => "hyperbolic triple-root capture",
    :K10 => "hyperbolic whirl in from infinity",
    :C_X1 => "parabolic double-horizon capture",
    :C_X2 => "hyperbolic double-horizon capture",
    :C_X3 => "parabolic triple-horizon capture",
    :C_X4 => "hyperbolic triple-horizon capture",
    :D1 => "parabolic scatter",
    :D2 => "hyperbolic scatter",
    :D_H1 => "parabolic horizon-root scatter",
    :D_H2 => "hyperbolic horizon-root scatter",
    :D_X1 => "parabolic extremal scatter",
    :D_X2 => "hyperbolic extremal scatter",
)

# Catalogue entries are keyed by the primary case ID (A1-A2, K1-K11, B1-B9, C1-C12, D1-D2,
# N1-N6) or, for the horizon-root (H) and extremal (X) members outside that numbering, by
# their tier name (A_H1: class letter, tier, number; displayed "A-H1", see
# `kerr_geo_case_name`). `broad_class` gives the class; example/data/catalogue_registry.tsv
# uses the same keys (display names).
Set(keys(KG56_MOTION_LABELS)) == Set(keys(KG56_REGISTRY)) ||
    error("KG56_MOTION_LABELS must cover exactly the 56 catalogue cases")

function _kg56_apex_constants(a, p, e, x)
    constants = kerr_geo_constants_of_motion(a, p, e, x)
    return (E=constants["E"], Lz=constants["Lz"], Q=constants["Q"])
end

function _kg56_critical_constants(a, energy, radius, branch)
    delta = radius^2 - 2radius + a^2
    root_term = sqrt(radius * (1 + (energy^2 - 1) * radius))
    momentum = delta * (energy * radius + branch * root_term) / (radius - 1)
    shifted = (energy * radius^2 - momentum) / a
    lz = a * energy + shifted
    k = momentum^2 / delta - radius^2
    return (E=energy, Lz=lz, Q=k - shifted^2)
end

function _kg56_example(id, broad, member_id, a, constants;
        member_slot, member_index=1, kwargs=(;), expected_case_ids,
        member_match_ids=(member_id,), view=:complete, sampling_window=nothing)
    return (
        case_id=Symbol(id),
        name=kerr_geo_case_name(Symbol(id)),
        broad_class=Symbol(broad),
        a=Float64(a),
        E=Float64(constants.E),
        Lz=Float64(constants.Lz),
        Q=Float64(constants.Q),
        member_slot=Symbol(member_slot),
        member_index=Int(member_index),
        kwargs=kwargs,
        expected_case_ids=Tuple(expected_case_ids),
        member_match_ids=Tuple(member_match_ids),
        view=Symbol(view),
        sampling_window=sampling_window,
    )
end

function _kg56_build_examples()
    examples = Any[]
    add(args...; kwargs...) = push!(examples, _kg56_example(args...; kwargs...))

    a_stable = 0.9
    b1 = _kg56_apex_constants(a_stable, 10.0, 0.5, 0.8)
    b2 = _kg56_apex_constants(a_stable, 8.0, 0.0, 0.8)
    isso = kerr_geo_isso(a_stable, 0.8)
    b3 = _kg56_apex_constants(a_stable, isso, 0.0, 0.8)
    ff17 = (
        E=0.9171300256198305,
        Lz=2.2591913519439517,
        Q=2.898994491984013,
    )
    parabolic_critical = (E=1.0, Lz=-0.7, Q=16.0)
    hyperbolic_critical = _kg56_critical_constants(0.7, 1.2, 3.0, 1.0)

    a_hc = 0.9
    rplus_hc = 1 + sqrt(1 - a_hc^2)
    hc_libration_E = 0.94
    hc_stable_E = 0.9331128235426359

    add(:A1, :stable, :A1, a_stable, b1;
        member_slot=:stable, expected_case_ids=(:A1, :B1))
    add(:A2, :stable, :A2, a_stable, b2;
        member_slot=:stable, expected_case_ids=(:A2, :B2))
    add(:A_H1, :stable, :A_H1, a_hc,
        (E=hc_libration_E, Lz=2rplus_hc * hc_libration_E / a_hc, Q=0.0);
        member_slot=:stable, expected_case_ids=(:A_H1,))
    add(:A_H2, :stable, :A_H2, a_hc,
        (E=hc_stable_E, Lz=2rplus_hc * hc_stable_E / a_hc, Q=0.0);
        member_slot=:stable, expected_case_ids=(:A_H2,))
    add(:A_X1, :stable, :A_X1, 1.0,
        (E=0.8, Lz=1.6, Q=1.0);
        member_slot=:stable, expected_case_ids=(:A_X1,))
    add(:A_X2, :stable, :A_X2, 1.0,
        (E=0.8, Lz=1.6, Q=0.8^4 / (1 - 0.8^2));
        member_slot=:stable, expected_case_ids=(:A_X2,))
    add(:K1, :critical, :K1, a_stable, b3;
        member_slot=:critical, member_index=1, expected_case_ids=(:K1, :K2))
    add(:K2, :critical, :K2, a_stable, b3;
        member_slot=:critical, member_index=2, expected_case_ids=(:K1, :K2))
    add(:K3, :critical, :K3, 0.7, ff17;
        member_slot=:critical, member_index=1, expected_case_ids=(:K3, :K4, :K5))
    add(:K4, :critical, :K4, 0.7, ff17;
        member_slot=:critical, member_index=2, expected_case_ids=(:K3, :K4, :K5))
    add(:K5, :critical, :K5, 0.7, ff17;
        member_slot=:critical, member_index=3, expected_case_ids=(:K3, :K4, :K5))
    add(:K6, :critical, :K6, 0.7, parabolic_critical;
        member_slot=:critical, member_index=1, expected_case_ids=(:K6, :K7, :K8))
    add(:K7, :critical, :K7, 0.7, parabolic_critical;
        member_slot=:critical, member_index=2, expected_case_ids=(:K6, :K7, :K8))
    add(:K8, :critical, :K8, 0.7, parabolic_critical;
        member_slot=:critical, member_index=3, expected_case_ids=(:K6, :K7, :K8))
    add(:K9, :critical, :K9, 0.7, hyperbolic_critical;
        member_slot=:critical, member_index=1, expected_case_ids=(:K9, :K10, :K11))
    add(:K10, :critical, :K10, 0.7, hyperbolic_critical;
        member_slot=:critical, member_index=2, expected_case_ids=(:K9, :K10, :K11))
    add(:K11, :critical, :K11, 0.7, hyperbolic_critical;
        member_slot=:critical, member_index=3, expected_case_ids=(:K9, :K10, :K11))
    add(:B1, :plunge, :B1, a_stable, b1;
        member_slot=:plunge, expected_case_ids=(:A1, :B1))
    add(:B2, :plunge, :B2, a_stable, b2;
        member_slot=:plunge, expected_case_ids=(:A2, :B2))
    add(:B3, :plunge, :B3, 0.9,
        (E=0.8105062006741774, Lz=1.1365738765215694,
         Q=0.0313190420592252);
        member_slot=:plunge, expected_case_ids=(:B3,))
    add(:B4, :plunge, :B4, 0.9, (E=0.94, Lz=0.1, Q=12.0);
        member_slot=:plunge, expected_case_ids=(:B4,))
    add(:B5, :plunge, :B5, 0.7, (E=1.0, Lz=4.0, Q=3.0);
        member_slot=:plunge, expected_case_ids=(:B5, :D1))
    add(:B6, :plunge, :B6, 0.5, (E=1.1, Lz=5.0, Q=1.0);
        member_slot=:plunge, expected_case_ids=(:B6, :D2))
    add(:B7, :plunge, :B7, 0.5,
        (E=0.5, Lz=0.20611404019925345, Q=7.611364448469808e-5);
        member_slot=:plunge, expected_case_ids=(:B7,))
    add(:B8, :plunge, :B8, 0.5,
        (E=0.5, Lz=0.19009568057487902, Q=0.0002955693714467326);
        member_slot=:plunge, expected_case_ids=(:B8,))
    add(:B9, :plunge, :B9, 0.0, (E=0.5, Lz=0.0, Q=0.0);
        member_slot=:plunge, kwargs=(axis=:north,),
        expected_case_ids=(:B9,))
    add(:B_X1, :plunge, :B_X1, 1.0, (E=0.8, Lz=1.6, Q=0.1);
        member_slot=:plunge, expected_case_ids=(:B_X1,),
        view=:future_plunge_half, sampling_window=(0.0, 8.0))
    add(:B_X2, :plunge, :B_X2, 1.0, (E=0.8, Lz=1.6, Q=0.92);
        member_slot=:plunge, expected_case_ids=(:B_X2,),
        view=:future_plunge_half, sampling_window=(0.0, 8.0))
    add(:C1, :capture, :C1, 0.5, (E=1.0, Lz=1.0, Q=0.0);
        member_slot=:capture, expected_case_ids=(:C1,))
    add(:C2, :capture, :C2, 0.9,
        (E=1.0, Lz=1.250039120491465, Q=0.005571204976471838);
        member_slot=:capture, expected_case_ids=(:C2,))
    add(:C3, :capture, :C3, 0.5, (E=1.1, Lz=0.2, Q=0.0);
        member_slot=:capture, expected_case_ids=(:C3,))
    add(:C4, :capture, :C4, 0.9,
        (E=1.1017435938802256, Lz=1.2421098461995382,
         Q=0.0014181133633207713);
        member_slot=:capture, expected_case_ids=(:C4,))
    add(:C5, :capture, :C5, 0.9, (E=1.8, Lz=0.2, Q=-1.0);
        member_slot=:capture, expected_case_ids=(:C5,))
    add(:C6, :capture, :C6, 0.5,
        (E=1.0, Lz=0.45836363636363636, Q=6.806611570247915e-5);
        member_slot=:capture, expected_case_ids=(:C6,))
    add(:C7, :capture, :C7, 0.5,
        (E=1.0, Lz=0.4492630807223752, Q=1.8558744446053963e-5);
        member_slot=:capture, expected_case_ids=(:C7,))
    add(:C8, :capture, :C8, 0.0, (E=1.0, Lz=0.0, Q=0.0);
        member_slot=:capture, kwargs=(axis=:north,),
        expected_case_ids=(:C8,))
    add(:C9, :capture, :C9, 0.5,
        (E=1.1, Lz=0.5087997013554778, Q=6.655227534990921e-5);
        member_slot=:capture, expected_case_ids=(:C9,))
    add(:C10, :capture, :C10, 0.5,
        (E=1.1, Lz=0.6456608711354802, Q=0.00162515588137737);
        member_slot=:capture, expected_case_ids=(:C10,))
    add(:C11, :capture, :C11, 0.5,
        (E=3.0, Lz=-0.01355932537398763, Q=-1.261626024149881);
        member_slot=:capture, expected_case_ids=(:C11,))
    add(:C12, :capture, :C12, 0.0, (E=1.5, Lz=0.0, Q=0.0);
        member_slot=:capture, kwargs=(axis=:north,),
        expected_case_ids=(:C12,))
    add(:C_X1, :capture, :C_X1, 1.0,
        (E=1.0, Lz=2.0, Q=1.0);
        member_slot=:capture, expected_case_ids=(:C_X1,))
    add(:C_X2, :capture, :C_X2, 1.0,
        (E=1.2, Lz=2.4, Q=1.0);
        member_slot=:capture, expected_case_ids=(:C_X2,))
    add(:C_X3, :capture, :C_X3, 1.0,
        (E=1.0, Lz=2.0, Q=2.0);
        member_slot=:capture, expected_case_ids=(:C_X3,))
    add(:C_X4, :capture, :C_X4, 1.0,
        (E=1.2, Lz=2.4, Q=3.32);
        member_slot=:capture, expected_case_ids=(:C_X4,))
    add(:D1, :scatter, :D1, 0.7, (E=1.0, Lz=4.0, Q=3.0);
        member_slot=:scatter, expected_case_ids=(:B5, :D1))
    add(:D2, :scatter, :D2, 0.5, (E=1.1, Lz=5.0, Q=1.0);
        member_slot=:scatter, expected_case_ids=(:B6, :D2))
    add(:D_H1, :scatter, :D_H1, a_hc,
        (E=1.0, Lz=2rplus_hc / a_hc, Q=0.0);
        member_slot=:scatter,
        expected_case_ids=(:D_H1,))
    add(:D_H2, :scatter, :D_H2, a_hc,
        (E=1.2, Lz=2rplus_hc * 1.2 / a_hc, Q=1.0);
        member_slot=:scatter,
        expected_case_ids=(:D_H2,))
    add(:D_X1, :scatter, :D_X1, 1.0,
        (E=1.0, Lz=2.0, Q=3.0);
        member_slot=:scatter,
        expected_case_ids=(:D_X1,))
    add(:D_X2, :scatter, :D_X2, 1.0,
        (E=1.2, Lz=2.4, Q=4.0);
        member_slot=:scatter, expected_case_ids=(:D_X2,))

    trapped_cases = (
        (:N1, :NFD01, 0.9, -0.95, -3.35, 0.0),
        (:N2, :NFD02, 0.9, -0.9409087099329637,
         -3.1442752637877538, 0.0),
        (:N3, :NFD03, 0.99, -0.2, -0.5, 0.0),
        (:N4, :NFD04, 0.9, -0.8, -4.0, 1.0),
        (:N5, :NFD05, 0.9, -1.0, -4.0, 0.0),
        (:N6, :NFD06, 0.9, -1.2, -5.0, 1.0),
    )
    for (id, disposition, a, E, Lz, Q) in trapped_cases
        add(id, :trapped, disposition, a, (E=E, Lz=Lz, Q=Q);
            member_slot=:trapped,
            kwargs=(trapped_component=:incoming, disposition_id=disposition),
            expected_case_ids=(id,), member_match_ids=(id, disposition),
            view=:future_plunge_half)
    end

    return Tuple(examples)
end

const KG56_EXAMPLES = _kg56_build_examples()
const KG56_EXAMPLE_BY_ID = Dict(example.case_id => example for example in KG56_EXAMPLES)

function kg56_example(case_id)
    id = kerr_geo_case_symbol(case_id)     # "D-H1" -> :D_H1
    haskey(KG56_EXAMPLE_BY_ID, id) || error("Unknown catalogue ID $(id). Use a case ID " *
        "(A1-A2, K1-K11, B1-B9, C1-C12, D1-D2, N1-N6) or a tier name (A-H1, B-X1, D-H1, ...).")
    return KG56_EXAMPLE_BY_ID[id]
end

function _kg56_case_ids(family)
    haskey(family.Status, :case_ids) || return ()
    return Tuple(family.Status.case_ids)
end

# the member's case ID and, for Critical members, the IDs of its constants' family
_kg56_member_ids(member) =
    Tuple(unique((member.CaseId, get(member.Status, :family_member_case_ids, ())...)))

function _kg56_select_member(family, example)
    slot = getfield(family, kerr_geo_class(example.broad_class).slot)
    member = slot isa Tuple ?
        (length(slot) >= example.member_index ? slot[example.member_index] : nothing) : slot
    member === nothing && error("$(example.case_id): no $(example.member_slot) member")
    return member
end

function kg56_construct(case_id)
    example = kg56_example(case_id)
    family = kerr_geodesic(
        example.a,
        (example.E, example.Lz, example.Q);
        example.kwargs...,
    )
    member = _kg56_select_member(family, example)
    return (example=example, family=family, member=member)
end

_kg56_trajectory(member) = member.Trajectory

function _kg56_callable(trajectory, names)
    for name in names
        haskey(trajectory, name) && return trajectory[name]
    end
    return nothing
end

_kg56_domain(member) = Tuple(member.Domain.mino)
_kg56_endpoint_roles(member) = Tuple(member.Domain.endpoint_roles)
_kg56_endpoint_closed(member) = Tuple(member.Domain.endpoint_closed)

function _kg56_uses_ingoing_coordinates(member, view)
    view === :future_plunge_half && return true
    return last(_kg56_endpoint_roles(member)) in (:future_horizon, :future_horizon_root_asymptote)
end

_kg56_spin(member) = member.ConstantsOfMotion.a

# on the spin axis φ is a gauge choice (it stands in for ψ)
_kg56_axis_display_gauge(member) = get(member.Roots, :polar, (;)) isa NamedTuple &&
    get(get(member.Roots, :polar, (;)), :sector, nothing) === :axis_constant

function _kg56_sampling_window(member)
    lower, upper = _kg56_domain(member)
    if isfinite(lower) && isfinite(upper)
        width = upper - lower
        width > 0 || return (lower - 4.0, upper + 4.0)
        margin = max(1.0e-7, 2.0e-3 * width)
        return (lower + margin, upper - margin)
    elseif !isfinite(lower) && !isfinite(upper)
        return (-8.0, 8.0)
    elseif !isfinite(lower)
        return (upper - 8.0, upper - 1.0e-5)
    else
        return (lower + 1.0e-5, lower + 8.0)
    end
end

function _kg56_lambda_grid(lower, upper, sample_count; cluster_right=false)
    sample_count >= 12 || error("sample_count must be at least 12")
    cluster_right || return collect(range(lower, upper; length=sample_count))
    sample_count < 24 && return collect(range(lower, upper; length=sample_count))

    # Near a repeated-root asymptote r(lambda) reaches the root to within rounding
    # exponentially fast, so a moderate global grid covers the branch; the remaining
    # probes go to its last quarter, where v(lambda) resolves the approach to the
    # future horizon.
    coarse_count = min(401, max(41, cld(sample_count, 10)))
    fine_count = sample_count - coarse_count + 1
    split = lower + 0.75 * (upper - lower)
    coarse = collect(range(lower, split; length=coarse_count))
    unit = collect(range(0.0, 1.0; length=fine_count))
    fine = split .+ (upper - split) .* (1 .- (1 .- unit) .^ 2)
    return vcat(coarse, fine[2:end])
end

function _kg56_past_repeated_asymptote(member)
    roles = _kg56_endpoint_roles(member)
    isempty(roles) && return false
    return occursin("repeated_root_asymptote", String(first(roles)))
end

function _kg56_future_repeated_asymptote(member)
    roles = _kg56_endpoint_roles(member)
    isempty(roles) && return false
    return occursin("repeated_root_asymptote", String(last(roles)))
end

function kg56_sample(member; sample_count=321, radial_cutoff=30.0,
        sampling_window=nothing, view=:complete)
    trajectory = _kg56_trajectory(member)
    rfun = _kg56_callable(trajectory, (:r,))
    thetafun = _kg56_callable(trajectory, (:theta,))
    zfun = _kg56_callable(trajectory, (:z,))
    phifun = _kg56_callable(trajectory, (:phi,))
    tfun = _kg56_callable(trajectory, (:t,))
    vfun = _kg56_callable(trajectory, (:v,))
    psifun = _kg56_callable(trajectory, (:psi,))
    rfun === nothing && error("Trajectory has no r(lambda) callable")
    (thetafun === nothing && zfun === nothing) &&
        error("Trajectory has neither theta(lambda) nor z(lambda)")
    phifun === nothing && error("Trajectory has no phi(lambda) callable")

    lower, upper = sampling_window === nothing ?
        _kg56_sampling_window(member) : sampling_window
    view === :future_plunge_half && (lower = max(lower, 0.0))
    use_ingoing = _kg56_uses_ingoing_coordinates(member, view)
    domain_upper = last(_kg56_domain(member))
    endpoint_closed = _kg56_endpoint_closed(member)
    if use_ingoing && isfinite(domain_upper) && endpoint_closed[2]
        upper = domain_upper
    end
    lower < upper || error("Sampling window must have positive width")

    parameter_name = use_ingoing ? :v : :t
    azimuth_name = use_ingoing ? :psi : :phi
    parameter_fun = use_ingoing ? vfun : tfun
    azimuth_fun = use_ingoing ? psifun : phifun
    if use_ingoing && azimuth_fun === nothing
        spin = _kg56_spin(member)
        (isfinite(spin) && abs(spin) <= 1.0e-14) ||
            _kg56_axis_display_gauge(member) ||
            error("Kerr horizon trajectory has no psi(lambda) callable")
        azimuth_name = abs(spin) <= 1.0e-14 ? :phi_equals_psi : :axis_phi_gauge
        azimuth_fun = phifun
    end
    parameter_fun === nothing &&
        error("Trajectory has no $(parameter_name)(lambda) callable")
    azimuth_fun === nothing &&
        error("Trajectory has no $(azimuth_name)(lambda) callable")

    lambda = Float64[]
    radius = Float64[]
    theta = Float64[]
    azimuth = Float64[]
    frame_parameter = Float64[]
    lambda_grid = _kg56_lambda_grid(
        lower,
        upper,
        sample_count;
        cluster_right=use_ingoing,
    )
    for value in lambda_grid
        try
            r = Float64(rfun(value))
            th = thetafun === nothing ?
                acos(clamp(Float64(zfun(value)), -1.0, 1.0)) :
                Float64(thetafun(value))
            az = Float64(azimuth_fun(value))
            parameter = Float64(parameter_fun(value))
            all(isfinite, (r, th, az, parameter)) || continue
            r > 0 || continue
            r <= radial_cutoff || continue
            push!(lambda, value)
            push!(radius, r)
            push!(theta, th)
            push!(azimuth, az)
            push!(frame_parameter, parameter)
        catch
        end
    end
    length(lambda) >= 12 || error(
        "Only $(length(lambda)) finite trajectory samples remained inside r <= $(radial_cutoff)")
    trimmed_past_asymptote_probes = 0
    trimmed_future_asymptote_probes = 0
    past_asymptote_radial_resolution = 0.0
    future_asymptote_radial_resolution = 0.0
    if _kg56_past_repeated_asymptote(member)
        radial_resolution = 1.0e-7 * max(1.0, abs(first(radius)))
        past_asymptote_radial_resolution = radial_resolution
        resolved = findfirst(
            index -> abs(radius[index] - first(radius)) > radial_resolution,
            eachindex(radius),
        )
        resolved === nothing && error(
            "The repeated-root departure is unresolved in the display window")
        trimmed_past_asymptote_probes = resolved - 1
        lambda = lambda[resolved:end]
        radius = radius[resolved:end]
        theta = theta[resolved:end]
        azimuth = azimuth[resolved:end]
        frame_parameter = frame_parameter[resolved:end]
    end
    if _kg56_future_repeated_asymptote(member)
        radial_resolution = 1.0e-7 * max(1.0, abs(last(radius)))
        future_asymptote_radial_resolution = radial_resolution
        resolved = findlast(
            index -> abs(radius[index] - last(radius)) > radial_resolution,
            eachindex(radius),
        )
        resolved === nothing && error(
            "The future repeated-root departure is unresolved in the display window")
        trimmed_future_asymptote_probes = length(radius) - resolved
        lambda = lambda[1:resolved]
        radius = radius[1:resolved]
        theta = theta[1:resolved]
        azimuth = azimuth[1:resolved]
        frame_parameter = frame_parameter[1:resolved]
    end
    length(lambda) >= 12 || error(
        "Only $(length(lambda)) resolved physical-time samples remain")
    all(diff(frame_parameter) .> 0) || error(
        "$(parameter_name)(lambda) is not strictly increasing on the display branch")
    x = radius .* sin.(theta) .* cos.(azimuth)
    y = radius .* sin.(theta) .* sin.(azimuth)
    z = radius .* cos.(theta)
    return (
        lambda=lambda,
        r=radius,
        theta=theta,
        azimuth=azimuth,
        frame_parameter=frame_parameter,
        frame_parameter_name=parameter_name,
        azimuth_name=azimuth_name,
        ingoing_coordinates=use_ingoing,
        trimmed_past_asymptote_probes=trimmed_past_asymptote_probes,
        trimmed_future_asymptote_probes=trimmed_future_asymptote_probes,
        past_asymptote_radial_resolution=past_asymptote_radial_resolution,
        future_asymptote_radial_resolution=
            future_asymptote_radial_resolution,
        x=x,
        y=y,
        z=z,
        requested_window=(lower, upper),
        radial_cutoff=radial_cutoff,
    )
end

function _kg56_radial_residual(member, lambda)
    try
        return abs(Float64(member.Residuals.radial(lambda)))
    catch
        return NaN
    end
end

function kg56_sample_catalog(; sample_count=161, radial_cutoff=30.0)
    Set(example.case_id for example in KG56_EXAMPLES) == Set(keys(KG56_REGISTRY)) &&
        length(KG56_EXAMPLES) == 56 ||
        error("KG56_EXAMPLES must hold exactly the 56 catalogue cases")

    rows = Any[]
    for example in KG56_EXAMPLES
        try
            built = kg56_construct(example.case_id)
            case_ids = _kg56_case_ids(built.family)
            case_pass = Set(case_ids) == Set(example.expected_case_ids)
            member_ids = _kg56_member_ids(built.member)
            member_pass = !isempty(intersect(
                Set(member_ids), Set(example.member_match_ids)))
            sample = kg56_sample(
                built.member;
                sample_count=sample_count,
                radial_cutoff=radial_cutoff,
                sampling_window=example.sampling_window,
                view=example.view,
            )
            midpoint = sample.lambda[cld(length(sample.lambda), 2)]
            radial_residual = _kg56_radial_residual(built.member, midpoint)
            pass = case_pass && member_pass && length(sample.lambda) >= 12
            push!(rows, (
                case_id=example.case_id,
                broad_class=example.broad_class,
                a=example.a,
                E=example.E,
                Lz=example.Lz,
                Q=example.Q,
                case_ids=case_ids,
                member_ids=member_ids,
                requested_sample_count=sample_count,
                sample_count=length(sample.lambda),
                trimmed_past_asymptote_probes=
                    sample.trimmed_past_asymptote_probes,
                trimmed_future_asymptote_probes=
                    sample.trimmed_future_asymptote_probes,
                min_r=minimum(sample.r),
                max_r=maximum(sample.r),
                frame_parameter=sample.frame_parameter_name,
                azimuth_parameter=sample.azimuth_name,
                frame_parameter_start=first(sample.frame_parameter),
                frame_parameter_end=last(sample.frame_parameter),
                endpoint_roles=_kg56_endpoint_roles(built.member),
                radial_residual=radial_residual,
                pass=pass,
                error="",
                sample=sample,
            ))
            built = nothing
            GC.gc(false)
        catch error
            push!(rows, (
                case_id=example.case_id,
                broad_class=example.broad_class,
                a=example.a,
                E=example.E,
                Lz=example.Lz,
                Q=example.Q,
                case_ids=(),
                member_ids=(),
                requested_sample_count=sample_count,
                sample_count=0,
                trimmed_past_asymptote_probes=0,
                trimmed_future_asymptote_probes=0,
                min_r=NaN,
                max_r=NaN,
                frame_parameter=:missing,
                azimuth_parameter=:missing,
                frame_parameter_start=NaN,
                frame_parameter_end=NaN,
                endpoint_roles=(),
                radial_residual=NaN,
                pass=false,
                error=sprint(showerror, error),
                sample=nothing,
            ))
            GC.gc(false)
        end
    end
    return rows
end

function _kg56_escape(value)
    text = string(value)
    text = replace(text, "&" => "&amp;")
    text = replace(text, "<" => "&lt;")
    text = replace(text, ">" => "&gt;")
    return replace(text, "\"" => "&quot;")
end

_kg56_number(value) = isfinite(value) ? @sprintf("%.8g", value) : "n/a"

function kg56_catalog_table()
    rows = String[]
    for example in KG56_EXAMPLES
        registry = KG56_REGISTRY[example.case_id]
        push!(rows, """
        <tr>
          <td><b>$(example.name)</b></td>
          <td>$(kerr_geo_class(example.broad_class).name)</td>
          <td>$(_kg56_number(example.a))</td>
          <td>$(_kg56_number(example.E))</td>
          <td>$(_kg56_number(example.Lz))</td>
          <td>$(_kg56_number(example.Q))</td>
          <td><code>$(_kg56_escape(registry.root_structure))</code></td>
          <td><code>$(_kg56_escape(registry.allowed_interval))</code></td>
          <td><code>$(_kg56_escape(registry.formula_family))</code></td>
        </tr>
        """)
    end
    return Base.HTML("""
    <div style="max-height:680px;overflow:auto;border:1px solid #d0d7de">
    <table style="border-collapse:collapse;width:100%;font:13px system-ui">
      <thead style="position:sticky;top:0;background:#f6f8fa">
        <tr>
          <th>ID</th><th>Class</th>
          <th>a</th><th>E</th><th>Lz</th><th>Q</th>
          <th>Roots</th><th>Allowed interval</th><th>Formula</th>
        </tr>
      </thead>
      <tbody>$(join(rows))</tbody>
    </table>
    </div>
    <style>
      table th, table td { padding:5px 7px; border-bottom:1px solid #d8dee4;
                           text-align:left; vertical-align:top; }
    </style>
    """)
end

function kg56_sample_table(catalog)
    rows = String[]
    for row in catalog
        color = row.pass ? "#1a7f37" : "#cf222e"
        status = row.pass ? "PASS" : "FAIL"
        residual = isnan(row.radial_residual) ? "not exposed" :
            _kg56_number(row.radial_residual)
        detail = row.pass ? "" : "<br><small>$(_kg56_escape(row.error))</small>"
        push!(rows, """
        <tr>
          <td><b>$(kerr_geo_case_name(row.case_id))</b></td>
          <td style="color:$(color);font-weight:700">$(status)</td>
          <td><code>$(_kg56_escape(row.case_ids))</code></td>
          <td>$(row.sample_count) / $(row.requested_sample_count)</td>
          <td>$(row.trimmed_past_asymptote_probes) / $(row.trimmed_future_asymptote_probes)</td>
          <td>$(_kg56_number(row.min_r))</td>
          <td>$(_kg56_number(row.max_r))</td>
          <td><code>$(row.frame_parameter) / $(row.azimuth_parameter)</code></td>
          <td>$(residual)$(detail)</td>
        </tr>
        """)
    end
    passed = count(row -> row.pass, catalog)
    summary = passed == length(catalog) ? "All $(passed) orbits built and sampled." :
        "$(passed) of $(length(catalog)) orbits built and sampled."
    return Base.HTML("""
    <p><b>$(summary)</b></p>
    <div style="max-height:620px;overflow:auto;border:1px solid #d0d7de">
    <table style="border-collapse:collapse;width:100%;font:13px system-ui">
      <thead style="position:sticky;top:0;background:#f6f8fa">
        <tr><th>ID</th><th>Status</th><th>Case IDs</th>
        <th>Retained / requested</th><th>Trimmed past / future</th>
        <th>min r</th><th>max r</th>
        <th>Frames / azimuth</th><th>Radial residual</th></tr>
      </thead>
      <tbody>$(join(rows))</tbody>
    </table></div>
    <style>
      table th, table td { padding:5px 7px; border-bottom:1px solid #d8dee4;
                           text-align:left; vertical-align:top; }
    </style>
    """)
end

function _kg56_polyline(u, v, x0, y0, width, height;
        umin=minimum(u), umax=maximum(u), vmin=minimum(v), vmax=maximum(v))
    urange = max(umax - umin, eps(Float64))
    vrange = max(vmax - vmin, eps(Float64))
    points = String[]
    for (ui, vi) in zip(u, v)
        x = x0 + width * (ui - umin) / urange
        y = y0 + height * (1 - (vi - vmin) / vrange)
        push!(points, @sprintf("%.2f,%.2f", x, y))
    end
    return join(points, " ")
end

function kg56_orbit_plot(case_id; sample_count=481, radial_cutoff=30.0)
    built = kg56_construct(case_id)
    example = built.example
    registry = KG56_REGISTRY[example.case_id]
    sample = kg56_sample(
        built.member;
        sample_count=sample_count,
        radial_cutoff=radial_cutoff,
        sampling_window=example.sampling_window,
        view=example.view,
    )
    rplus = 1 + sqrt(max(0.0, 1 - example.a^2))
    spatial_limit = max(
        1.15rplus,
        maximum(abs, vcat(sample.x, sample.y, sample.z)),
    )
    spatial_limit = min(spatial_limit, radial_cutoff)
    xy = _kg56_polyline(
        sample.x, sample.y, 55, 62, 245, 225;
        umin=-spatial_limit, umax=spatial_limit,
        vmin=-spatial_limit, vmax=spatial_limit,
    )
    xz = _kg56_polyline(
        sample.x, sample.z, 365, 62, 245, 225;
        umin=-spatial_limit, umax=spatial_limit,
        vmin=-spatial_limit, vmax=spatial_limit,
    )
    rr = _kg56_polyline(
        sample.lambda, sample.r, 675, 62, 245, 225;
        vmin=minimum(sample.r), vmax=maximum(sample.r),
    )
    horizon_radius = 245 * rplus / (2spatial_limit)
    roles = join(replace.(string.(built.member.Domain.endpoint_roles), "_" => " "), " → ")
    title = "$(example.name): $(KG56_MOTION_LABELS[example.case_id]), $(roles)"
    parameters = "a=$(_kg56_number(example.a)), E=$(_kg56_number(example.E)), " *
        "Lz=$(_kg56_number(example.Lz)), Q=$(_kg56_number(example.Q))"
    svg = """
    <div style="font:14px system-ui;max-width:980px">
      <h3 style="margin:8px 0 2px">$(_kg56_escape(title))</h3>
      <div style="color:#57606a;margin-bottom:7px">$(_kg56_escape(parameters))</div>
      <svg viewBox="0 0 960 330" width="100%" role="img"
           aria-label="$(_kg56_escape(title))">
        <rect x="0" y="0" width="960" height="330" fill="#ffffff"/>
        <g stroke="#d0d7de" stroke-width="1" fill="none">
          <rect x="55" y="62" width="245" height="225"/>
          <rect x="365" y="62" width="245" height="225"/>
          <rect x="675" y="62" width="245" height="225"/>
          <line x1="177.5" y1="62" x2="177.5" y2="287"/>
          <line x1="55" y1="174.5" x2="300" y2="174.5"/>
          <line x1="487.5" y1="62" x2="487.5" y2="287"/>
          <line x1="365" y1="174.5" x2="610" y2="174.5"/>
        </g>
        <g fill="#24292f" font-family="system-ui" font-size="13">
          <text x="55" y="48">x-y projection</text>
          <text x="365" y="48">x-z projection</text>
          <text x="675" y="48">radial motion r(lambda)</text>
          <text x="55" y="309">spatial range: +/-$(_kg56_number(spatial_limit)) M</text>
          <text x="675" y="309">lambda in [$(_kg56_number(first(sample.lambda))), $(_kg56_number(last(sample.lambda)))]</text>
        </g>
        <circle cx="177.5" cy="174.5" r="$(@sprintf("%.2f", horizon_radius))"
                fill="#111827" fill-opacity="0.72"/>
        <circle cx="487.5" cy="174.5" r="$(@sprintf("%.2f", horizon_radius))"
                fill="#111827" fill-opacity="0.72"/>
        <polyline points="$(xy)" fill="none" stroke="#0969da" stroke-width="2"/>
        <polyline points="$(xz)" fill="none" stroke="#8250df" stroke-width="2"/>
        <polyline points="$(rr)" fill="none" stroke="#cf222e" stroke-width="2"/>
      </svg>
      <div style="color:#57606a">
        Formula: <code>$(_kg56_escape(registry.formula_family))</code>;
        allowed interval: <code>$(_kg56_escape(registry.allowed_interval))</code>.
      </div>
    </div>
    """
    return Base.HTML(svg)
end

function kg56_case_summary(case_id)
    built = kg56_construct(case_id)
    example = built.example
    registry = KG56_REGISTRY[example.case_id]
    case_ids = _kg56_case_ids(built.family)
    member_ids = _kg56_member_ids(built.member)
    return Base.HTML("""
    <table style="border-collapse:collapse;font:14px system-ui">
      <tr><th>Case</th><td><b>$(example.name)</b> ($(kerr_geo_class(example.broad_class).name))</td></tr>
      <tr><th>Input</th><td><code>(a=$(example.a), E=$(example.E), Lz=$(example.Lz), Q=$(example.Q))</code></td></tr>
      <tr><th>Case IDs</th><td><code>$(_kg56_escape(case_ids))</code></td></tr>
      <tr><th>Selected member</th><td><code>$(example.member_slot)[$(example.member_index)]</code></td></tr>
      <tr><th>Member IDs</th><td><code>$(_kg56_escape(member_ids))</code></td></tr>
      <tr><th>Root structure</th><td><code>$(_kg56_escape(registry.root_structure))</code></td></tr>
      <tr><th>Allowed interval</th><td><code>$(_kg56_escape(registry.allowed_interval))</code></td></tr>
      <tr><th>Formula family</th><td><code>$(_kg56_escape(registry.formula_family))</code></td></tr>
    </table>
    <style>table th,table td{padding:6px 9px;border-bottom:1px solid #d8dee4;text-align:left}</style>
    """)
end

function kg56_custom_classification(a, E, Lz, Q; kwargs...)
    family = kerr_geodesic(a, (E, Lz, Q); kwargs...)
    case_ids = _kg56_case_ids(family)
    return (
        input=(a=a, E=E, Lz=Lz, Q=Q),
        case_names=Tuple(kerr_geo_case_name.(case_ids)),
        case_ids=case_ids,
        family=family,
    )
end

function kg56_gallery(catalog)
    cards = String[]
    for row in catalog
        motion_label = KG56_MOTION_LABELS[row.case_id]
        color = row.pass ? "#0969da" : "#cf222e"
        body = if row.pass
            sample = row.sample
            limit = max(1.0, maximum(abs, vcat(sample.x, sample.y)))
            points = _kg56_polyline(
                sample.x, sample.y, 8, 8, 204, 150;
                umin=-limit, umax=limit, vmin=-limit, vmax=limit,
            )
            """<svg viewBox="0 0 220 166" width="100%">
                 <rect width="220" height="166" fill="#fff"/>
                 <polyline points="$(points)" fill="none" stroke="$(color)" stroke-width="1.8"/>
               </svg>"""
        else
            """<div style="height:150px;padding:8px;color:#cf222e">$(_kg56_escape(row.error))</div>"""
        end
        push!(cards, """
        <div style="border:1px solid #d0d7de;background:#fff">
          <div style="padding:7px 9px;border-bottom:1px solid #d8dee4">
            <div><b>$(kerr_geo_case_name(row.case_id))</b> <span style="color:#57606a">$(kerr_geo_class(row.broad_class).name)</span></div>
            <div style="margin-top:2px;font-weight:600">$(_kg56_escape(motion_label))</div>
          </div>
          $(body)
          <div style="padding:5px 9px;color:#57606a;font-size:12px">
            a=$(_kg56_number(row.a)), E=$(_kg56_number(row.E)),
            Lz=$(_kg56_number(row.Lz)), Q=$(_kg56_number(row.Q))
          </div>
        </div>
        """)
    end
    return Base.HTML("""
    <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));
                gap:10px;font:13px system-ui">$(join(cards))</div>
    """)
end

function _kg56_linear_interpolate(grid, values, target)
    target <= first(grid) && return first(values)
    target >= last(grid) && return last(values)
    left = searchsortedlast(grid, target)
    right = left + 1
    weight = (target - grid[left]) / (grid[right] - grid[left])
    return muladd(weight, values[right] - values[left], values[left])
end

function _kg56_uniform_frame_sample(sample, point_limit)
    count = min(max(2, Int(point_limit)), length(sample.frame_parameter))
    parameter = collect(range(
        first(sample.frame_parameter),
        last(sample.frame_parameter);
        length=count,
    ))
    interpolate(values) = [
        _kg56_linear_interpolate(sample.frame_parameter, values, target)
        for target in parameter
    ]
    return (
        parameter=parameter,
        lambda=interpolate(sample.lambda),
        x=interpolate(sample.x),
        y=interpolate(sample.y),
        z=interpolate(sample.z),
    )
end

# Animated gallery in the look of the README animation (example/animations/showcase_all.gif):
# dark sky, horizon with a soft glow, faint ergosphere wireframe and spin axis, the whole path
# faint and a fading trail in the class palette, everything behind the horizon hidden. Each
# camera is fixed at the elevation/azimuth of views.tsv (orthographic, matplotlib's view_init
# convention). Julia projects the geometry once; the
# static layer is SVG (it shows even where output scripts do not run) and a short script
# draws the moving trail and particle on a canvas above it.

const KG56_VIEWS_PATH = joinpath(KG56_NOTEBOOK_DIR, "views.tsv")
const KG56_DEFAULT_VIEW = (elev=30.0, azim=-60.0)
const KG56_GALLERY_COUNT = Ref(0)
const KG56_TRAIL_BINS = 20
const KG56_HIDDEN = typemin(Int16)     # projected point behind the horizon

function _kg56_read_views(path=KG56_VIEWS_PATH)
    views = Dict{String,NamedTuple{(:elev, :azim, :description),
        Tuple{Float64,Float64,String}}}()
    for line in eachline(path)
        (isempty(line) || startswith(line, "#") || startswith(line, "name\t")) && continue
        name, elev, azim, description = split(line, '\t'; limit=4)
        views[name] = (
            elev=isempty(elev) ? KG56_DEFAULT_VIEW.elev : parse(Float64, elev),
            azim=isempty(azim) ? KG56_DEFAULT_VIEW.azim : parse(Float64, azim),
            description=String(description),
        )
    end
    return views
end

kg56_bounded(endpoint_roles) =
    !any(role -> occursin("infinity", String(role)), endpoint_roles)

_kg56_rgb(hex) = Tuple(parse(Int, hex[i:i+1]; base=16) / 255 for i in (2, 4, 6))

# matplotlib's LinearSegmentedColormap.from_list: the palette stops evenly spaced on [0, 1]
function _kg56_colormap(palette, x)
    stops = map(_kg56_rgb, collect(palette))
    s = clamp(x, 0.0, 1.0) * (length(stops) - 1)
    k = min(floor(Int, s), length(stops) - 2)
    w = s - k
    return Tuple((1 - w) * stops[k+1][c] + w * stops[k+2][c] for c in 1:3)
end

function _kg56_css(rgb, alpha=1.0)
    r, g, b = round.(Int, 255 .* rgb)
    return @sprintf("rgba(%d,%d,%d,%.3f)", r, g, b, alpha)
end

# Panel coordinates: the view spans 1000 x 1000 units and the frame limit lies 357 units
# from the centre (what matplotlib's 3D axes give in the README animation: 0.3568 of the figure
# width); stored in tenths of a unit as Int16.
_kg56_panel_unit(value) = round(Int16, 10 * clamp(value, -3000.0, 3000.0))

function _kg56_svg_path(u, v; visible=trues(length(u)))
    io = IOBuffer()
    for i in eachindex(u)
        visible[i] || continue
        pen = i > firstindex(u) && visible[i-1] ? 'L' : 'M'
        print(io, pen, round(Int, u[i]), ' ', round(Int, v[i]))
    end
    return String(take!(io))
end

# Projected geometry of one panel, as in the README animation: positions uniform in the
# frame parameter (t or v), r - r+ magnified for orbits that stay near the horizon, the
# same frame limit, and `points` positions uniform in a blend of that time and arc length
# (the particle moves uniformly through them).
function _kg56_panel_geometry(sample, a, bounded, elev, azim; points=1600)
    dense = _kg56_uniform_frame_sample(sample, 6000)
    X, Y, Z = dense.x, dense.y, dense.z
    rplus = 1 + sqrt(max(0.0, 1 - a^2))
    r = sqrt.(X .^ 2 .+ Y .^ 2 .+ Z .^ 2)
    magnify = 1.0
    if maximum(r) - rplus < 0.8rplus
        magnify = 1.2rplus / max(maximum(r) - rplus, 1.0e-9)
        stretch = (rplus .+ magnify .* (r .- rplus)) ./ max.(r, 1.0e-12)
        X, Y, Z = X .* stretch, Y .* stretch, Z .* stretch
    end
    extent = maximum(abs, vcat(X, Y, Z))
    limit = max(3.2rplus, bounded ? 0.9extent + 0.4 : min(0.62extent + 0.4, 16.0))

    frame = dense.parameter .- first(dense.parameter)
    arc = vcat(0.0, cumsum(sqrt.(diff(X) .^ 2 .+ diff(Y) .^ 2 .+ diff(Z) .^ 2)))
    blend = 0.5 .* frame ./ max(last(frame), 1.0e-300) .+ 0.5 .* arc ./ max(last(arc), 1.0e-300)
    grid = range(0.0, 1.0; length=points)
    along(values) = [_kg56_linear_interpolate(blend, values, u) for u in grid]
    x, y, z, f = along(X), along(Y), along(Z), along(frame)

    se, ce = sincosd(elev)
    sa, ca = sincosd(azim)
    scale = 356.8 / limit
    horizontal(x, y, z) = -sa * x + ca * y
    vertical(x, y, z) = -se * ca * x - se * sa * y + ce * z
    depth(x, y, z) = ce * ca * x + ce * sa * y + se * z
    panel_u(x, y, z) = 500 + scale * horizontal(x, y, z)
    panel_v(x, y, z) = 500 - scale * vertical(x, y, z)
    function hidden(x, y, z)
        d2 = horizontal(x, y, z)^2 + vertical(x, y, z)^2
        return d2 < rplus^2 && depth(x, y, z) < sqrt(max(rplus^2 - d2, 0.0))
    end

    u, v = panel_u.(x, y, z), panel_v.(x, y, z)
    visible = .!hidden.(x, y, z)
    xy = Vector{Int16}(undef, 2points)
    for i in 1:points
        xy[2i-1] = visible[i] ? _kg56_panel_unit(u[i]) : KG56_HIDDEN
        xy[2i] = _kg56_panel_unit(v[i])
    end

    # ergosphere r = 1 + sqrt(1 - a^2 cos^2 theta): 20 meridians and 4 parallels, as the
    # wireframe of the README animation
    ergo(theta) = 1 + sqrt(max(1 - a^2 * cos(theta)^2, 0.0))
    shell(phi, theta) = (ergo(theta) * sin(theta) * cos(phi),
        ergo(theta) * sin(theta) * sin(phi), ergo(theta) * cos(theta))
    curves = vcat(
        [[shell(2pi * k / 20, theta) for theta in range(0, pi; length=20)] for k in 0:19],
        [[shell(phi, 4pi * k / 19) for phi in range(0, 2pi; length=41)] for k in 1:4],
    )
    wireframe = join(_kg56_svg_path(panel_u.(first.(c), getindex.(c, 2), last.(c)),
        panel_v.(first.(c), getindex.(c, 2), last.(c))) for c in curves)
    axis = _kg56_svg_path(panel_u.(0.0, 0.0, [-1.35rplus, 1.35rplus]),
        panel_v.(0.0, 0.0, [-1.35rplus, 1.35rplus]))

    return (
        xy=xy,
        frame=Float32.(f),
        tail=0.25 * last(frame),
        path=_kg56_svg_path(u, v; visible=visible),
        wireframe=wireframe,
        axis=axis,
        horizon=scale * rplus,
        glow=2.6 * 490 / limit * rplus,    # the glow radius of the README animation
        magnify=magnify,
        limit=limit,
    )
end

function _kg56_star_symbol(id; count=150)
    # fixed pseudo-random sky; each panel shows it rotated or mirrored
    noise(k, seed) = mod(sin(12.9898k + 78.233seed) * 43758.5453, 1.0)
    stars = join(@sprintf("<circle cx=\"%.0f\" cy=\"%.0f\" r=\"%.1f\" fill-opacity=\"%.2f\"/>",
        1000noise(k, 1), 1000noise(k, 2), 0.9 + 2.2noise(k, 3)^2, 0.18 + 0.3noise(k, 4))
        for k in 1:count)
    return """<symbol id="$(id)" viewBox="0 0 1000 1000"><g fill="#ffffff">$(stars)</g></symbol>"""
end

function _kg56_glow_stops()
    # alpha(rho) = 0.2 * clip((2.6 R - rho) / (1.56 R), 0, 1)^1.6 over 0 <= rho <= 2.6 R
    return join(@sprintf("<stop offset=\"%.2f\" stop-color=\"#6d8cf0\" stop-opacity=\"%.4f\"/>",
        o, 0.2 * clamp(2.6 * (1 - o) / 1.56, 0.0, 1.0)^1.6) for o in 0.0:0.05:1.0)
end

function _kg56_panel_card(index, row, geometry, view, gid)
    class = kerr_geo_class(row.broad_class)
    name = kerr_geo_case_name(row.case_id)
    stars = @sprintf("rotate(%d 500 500)%s", 90 * (index % 4),
        (index ÷ 4) % 2 == 1 ? " translate(1000 0) scale(-1 1)" : "")
    path_color = _kg56_css(_kg56_colormap(class.palette, 0.6))
    magnified = geometry.magnify > 1 ?
        """<div style="font-size:9px;opacity:.6;margin-top:1px">r &minus; r&#8330; magnified
        &times;$(@sprintf("%.3g", geometry.magnify))</div>""" : ""
    stroke = "fill=\"none\" vector-effect=\"non-scaling-stroke\" stroke-linecap=\"round\" stroke-linejoin=\"round\""
    constants = @sprintf("a = %.4g&ensp; E = %.4g&ensp; L<sub>z</sub> = %.4g&ensp; Q = %.4g",
        row.a, row.E, row.Lz, row.Q)
    return """
    <div data-kg56="$(index)" style="position:relative;aspect-ratio:1/1;background:#05070d;
         overflow:hidden;min-width:0">
      <svg viewBox="0 0 1000 1000" style="position:absolute;inset:0;width:100%;height:100%"
           role="img" aria-label="Animated trajectory for $(name): $(_kg56_escape(KG56_MOTION_LABELS[row.case_id]))">
        <use href="#$(gid)-stars" xlink:href="#$(gid)-stars" width="1000" height="1000"
             transform="$(stars)"/>
        <circle cx="500" cy="500" r="$(@sprintf("%.1f", geometry.glow))" fill="url(#$(gid)-glow)"/>
        <circle class="kg56-horizon" cx="500" cy="500" r="$(@sprintf("%.1f", geometry.horizon))"
                fill="url(#$(gid)-horizon)"/>
        <path d="$(geometry.wireframe)" stroke="#6f7fb0" stroke-opacity="0.22" stroke-width="0.5" $(stroke)/>
        <path d="$(geometry.axis)" stroke="#8fa3d9" stroke-opacity="0.4" stroke-width="0.7" $(stroke)/>
        <path class="kg56-path" d="$(geometry.path)" stroke="$(path_color)" stroke-opacity="0.34"
              stroke-width="0.8" $(stroke)/>
      </svg>
      <canvas style="position:absolute;inset:0;width:100%;height:100%"></canvas>
      <div style="position:absolute;left:4.5%;top:3.5%;right:4.5%;color:#dfe6ff;line-height:1.3">
        <div><b style="font-size:15px">$(name)</b>
          <span style="font-size:11px;opacity:.9">&ensp;$(_kg56_escape(KG56_MOTION_LABELS[row.case_id]))</span></div>
        <div style="font-size:10px;opacity:.85"><span style="color:$(class.color)">$(class.name)</span>
          &middot; $(_kg56_escape(view.description))</div>
        $(magnified)
      </div>
      <div style="position:absolute;left:4.5%;right:4.5%;bottom:3%;display:flex;
           justify-content:space-between;color:#dfe6ff;opacity:.7;font-size:9px">
        <span>$(constants)</span><span class="kg56-clock">$(row.sample.frame_parameter_name) = 0.0 M</span>
      </div>
    </div>
    """
end

function _kg56_panel_data(row, geometry)
    class = kerr_geo_class(row.broad_class)
    trail = join(("\"" * _kg56_css(_kg56_colormap(class.palette, 0.25 + 0.75w), w^1.6) * "\""
        for w in ((0:KG56_TRAIL_BINS-1) .+ 0.5) ./ KG56_TRAIL_BINS), ",")
    glow = _kg56_css(_kg56_rgb(class.color), 0.45)
    return """{"xy":"$(base64encode(geometry.xy))","f":"$(base64encode(geometry.frame))",""" *
        """"tail":$(geometry.tail),"frame":"$(row.sample.frame_parameter_name)",""" *
        """"trail":[$(trail)],"glow":"$(glow)"}"""
end

const KG56_GALLERY_SCRIPT = raw"""
(() => {
  const DATA = __KG56_DATA__, SECONDS = __KG56_SECONDS__;
  const me = document.currentScript;
  const root = (me && me.closest('[data-kg56-gallery]')) || document.getElementById(__KG56_ID__);
  if (!root) return;
  const bytes = b => Uint8Array.from(atob(b), c => c.charCodeAt(0)).buffer;
  const HIDDEN = -32768, BINS = DATA.length ? DATA[0].trail.length : 0;
  const cards = DATA.map((d, i) => {
    const el = root.querySelector('[data-kg56="' + i + '"]');
    const canvas = el && el.querySelector('canvas');
    if (!canvas) return null;
    return {d, canvas, ctx: canvas.getContext('2d'), clock: el.querySelector('.kg56-clock'),
            xy: new Int16Array(bytes(d.xy)), f: new Float32Array(bytes(d.f)), shown: true, text: ''};
  }).filter(Boolean);
  if ('IntersectionObserver' in window) {
    const seen = new IntersectionObserver(entries => entries.forEach(e => {
      e.target.kg56.shown = e.isIntersecting; }));
    cards.forEach(c => { c.canvas.kg56 = c; c.shown = false; seen.observe(c.canvas); });
  }

  function draw(c, phase) {
    const {canvas, ctx, xy, f, d} = c;
    const width = canvas.clientWidth;
    if (!width) return;
    const px = Math.round(width * (window.devicePixelRatio || 1));
    if (canvas.width !== px) { canvas.width = px; canvas.height = px; }
    const unit = 1000 / width, grow = Math.max(1, width / 320);   // one CSS pixel in panel units
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, px, px);
    ctx.setTransform(px / 1000, 0, 0, px / 1000, 0, 0);

    const n = f.length, k = phase * (n - 1), j = Math.min(Math.floor(k), n - 2), w1 = k - j;
    const X = i => xy[2 * i] / 10, Y = i => xy[2 * i + 1] / 10, hid = i => xy[2 * i] === HIDDEN;
    const now = f[j] + w1 * (f[j + 1] - f[j]), start = now - d.tail;
    let lo = 0, hi = j;                                   // first point of the trail
    while (lo < hi) { const m = (lo + hi) >> 1; if (f[m] < start) lo = m + 1; else hi = m; }
    const bins = Array.from({length: BINS}, () => new Path2D());
    for (let i = Math.max(lo - 1, 0); i <= j; i++) {
      if (hid(i) || hid(i + 1)) continue;
      const head = i === j, ax = X(i), ay = Y(i);
      const bx = head ? ax + w1 * (X(i + 1) - ax) : X(i + 1);
      const by = head ? ay + w1 * (Y(i + 1) - ay) : Y(i + 1);
      const mid = 0.5 * (f[i] + (head ? now : f[i + 1]));
      const w = Math.min(Math.max((mid - start) / Math.max(d.tail, 1e-300), 0), 1);
      const b = Math.min(BINS - 1, Math.floor(w * BINS));
      bins[b].moveTo(ax, ay); bins[b].lineTo(bx, by);
    }
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    for (let b = 0; b < BINS; b++) {
      ctx.strokeStyle = d.trail[b];
      ctx.lineWidth = (0.45 + 2.1 * (b + 0.5) / BINS) * grow * unit;
      ctx.stroke(bins[b]);
    }
    if (!hid(j) && !hid(j + 1)) {
      const x = X(j) + w1 * (X(j + 1) - X(j)), y = Y(j) + w1 * (Y(j + 1) - Y(j));
      const r = 7 * grow * unit, glow = ctx.createRadialGradient(x, y, 0, x, y, r);
      glow.addColorStop(0, d.glow); glow.addColorStop(1, 'rgba(0,0,0,0)');
      ctx.fillStyle = glow; ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = '#ffffff'; ctx.beginPath(); ctx.arc(x, y, 1.6 * grow * unit, 0, 2 * Math.PI);
      ctx.fill();
    }
    const text = d.frame + ' = ' + now.toFixed(1) + ' M';
    if (c.clock && text !== c.text) { c.clock.textContent = text; c.text = text; }
  }

  let attached = false, last = -Infinity;
  const t0 = performance.now();
  function frame(time) {
    if (root.isConnected) attached = true; else if (attached) return;   // output cleared
    requestAnimationFrame(frame);
    if (time - last < 1000 / 30 - 2) return;                           // 30 fps
    last = time;
    const phase = (((time - t0) / 1000) % SECONDS) / SECONDS;
    for (const c of cards) if (c.shown) draw(c, phase);
  }
  requestAnimationFrame(frame);
})();
"""

function kg56_animated_gallery(catalog; points=1600, seconds=8.0,
        views=_kg56_read_views())
    gid = "kg56-gallery-$(KG56_GALLERY_COUNT[] += 1)"
    cards = String[]
    data = String[]
    for row in catalog
        name = kerr_geo_case_name(row.case_id)
        if !(row.pass && row.sample !== nothing)
            push!(cards, """<div style="aspect-ratio:1/1;background:#05070d;color:#ff7b72;
                padding:10px;font-size:12px"><b>$(name)</b><br>$(_kg56_escape(row.error))</div>""")
            continue
        end
        view = get(views, name,
            merge(KG56_DEFAULT_VIEW, (description=KG56_MOTION_LABELS[row.case_id],)))
        geometry = _kg56_panel_geometry(row.sample, row.a, kg56_bounded(row.endpoint_roles),
            view.elev, view.azim; points=points)
        push!(cards, _kg56_panel_card(length(data), row, geometry, view, gid))
        push!(data, _kg56_panel_data(row, geometry))
    end
    legend = join(("""<span style="color:$(c.color)">&#9679; $(c.name)</span>"""
        for c in KERR_GEO_CLASSES), "&emsp;")
    script = replace(KG56_GALLERY_SCRIPT, "__KG56_ID__" => "\"$(gid)\"",
        "__KG56_DATA__" => "[$(join(data, ","))]", "__KG56_SECONDS__" => string(Float64(seconds)))
    return Base.HTML("""
    <div id="$(gid)" data-kg56-gallery style="background:#05070d;color:#c9d3f2;
         font:12px system-ui,-apple-system,sans-serif;padding:12px">
      <svg width="0" height="0" style="position:absolute" aria-hidden="true">
        <defs>
          $(_kg56_star_symbol("$(gid)-stars"))
          <radialGradient id="$(gid)-glow">$(_kg56_glow_stops())</radialGradient>
          <radialGradient id="$(gid)-horizon" cx="0.5" cy="0.5" r="0.5" fx="0.36" fy="0.32">
            <stop offset="0" stop-color="#232c50"/><stop offset="0.55" stop-color="#121832"/>
            <stop offset="1" stop-color="#090c19"/>
          </radialGradient>
        </defs>
      </svg>
      <div style="margin:0 2px 10px;line-height:1.45">
        Each panel replays one orbit with a fixed orthographic camera (elevation and azimuth
        from example/views.tsv), in the style of the README animation. The particle advances
        uniformly in a blend of the physical time and arc length; the trail covers the last
        quarter of that time. Horizon-ending branches use advanced time v and ingoing azimuth
        &psi;, all others t and &phi;. Orbits confined near the horizon are drawn with
        r &minus; r&#8330; magnified. Parts behind the horizon are hidden.<br>$(legend)
      </div>
      <div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(250px,1fr));
                  gap:3px">$(join(cards))</div>
      <script>$(script)</script>
    </div>
    """)
end
