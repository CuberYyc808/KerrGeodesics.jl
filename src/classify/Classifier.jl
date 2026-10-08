# All-component classifier: radial roots of (a, E, Lz, Q) -> allowed exterior intervals ->
# case IDs and KerrGeoRadialComponent records, plus polar admissibility and component selection.


const CLASSIFICATION_PIPELINE = (
    :canonical_constants,
    :metric_limit,
    :energy_regime,
    :radial_polynomial_degree,
    :roots_and_multiplicities,
    :allowed_radial_dynamical_components,
    :polar_admissibility,
    :future_direction,
    :case_assignment,
    :family_assembly,
    :optional_selection,
)

"""
    kerr_geo_classification_pipeline()

The stages of `kerr_geo_classify`, in order: canonical constants, metric limit, energy
regime, radial degree, roots and multiplicities, allowed components, polar admissibility,
future direction, case assignment, family assembly and the optional selection of one
component.
"""
kerr_geo_classification_pipeline() = CLASSIFICATION_PIPELINE

_root_close(x, y; atol=_root_atol(_float_type(x, y)), rtol=_root_rtol(_float_type(x, y))) =
    abs(x - y) <= atol + rtol * max(1.0, abs(x), abs(y))

# Whether a repeated root of multiplicity m sits at r, and where. Exterior (r ≥ r₊, or a
# cluster on the horizon itself when P(r₊) = 0): the roots bound the motion, so the constants
# decide (`_repeated_radial_zero`: within one ulp per component of the repeated-root
# manifold) and the root stays at r. Interior (r < r₊): the roots only enter the closed form
# of r(λ), evaluated at r ≥ r₊. Replacing the m estimates nearest r by one root x changes R
# there by the factor Π(r₊ − zᵢ)/(r₊ − x)^m, largest at r₊; with x the cluster mean the change
# is O(δ²/(r₊ − x)²) for a cluster of size δ. The cluster is merged when the change of
# R(r₊) is within the one-ulp input reach at r₊. Returns the
# radius and how the repeated root was read, or `nothing`:
#   :exact                  the constants have this repeated root (R, …, R^(m−1) vanish);
#   :within_input_ulp       nonexact repeated reading within the input reach.
function _accept_repeated(a, energy, lz, q, r, m, raw_roots, rplus, located)
    ordered = sort(raw_roots; by=z -> abs(z - r))
    k = min(m, length(ordered))
    nearest, rest = ordered[1:k], ordered[k+1:end]
    on_horizon = _horizon_root(a, energy, lz) && all(z -> abs(rplus - r) < abs(z - r), rest)
    if r >= rplus || on_horizon
        _repeated_radial_zero(a, energy, lz, q, r; multiplicity=m, located=located) || return nothing
        exact = _InputArithmetic.repeated_at(a, energy, lz, q, r, m)
        return (r, exact ? :exact : :within_input_ulp)
    end
    all(z -> real(z) < rplus, nearest) || return nothing
    # the located radius r (exact for an exact repeated root) or the mean, whichever is closer
    change(x) = abs(prod(rplus .- nearest) / (rplus - x)^m - 1)
    c = real(sum(nearest) / k)
    x = change(r) <= change(c) ? r : c
    edge = abs(first(_radial_root_evaluator(a, energy, lz, q)(rplus)))
    change(x) * edge <= _radial_input_reach(a, energy, lz, q, rplus, 0) || return nothing
    exact = _InputArithmetic.repeated_at(a, energy, lz, q, x, m)
    return (x, exact ? :exact : :within_input_ulp)
end

function _repeated_root_candidates(a, energy, lz, q, polynomial, raw_roots, rplus;
        atol=_root_atol(_float_type(a, energy, lz, q)), rtol=_root_rtol(_float_type(a, energy, lz, q)))
    candidates = NamedTuple[]
    diagnostics = NamedTuple[]
    derivative_polynomial = derivative(polynomial)
    for derivative_order in 1:max(0, degree(polynomial) - 1)
        for value in _polynomial_roots(coeffs(derivative_polynomial))
            imag_tolerance = atol + rtol * max(1.0, abs(real(value)))
            abs(imag(value)) <= imag_tolerance || continue
            radius = float(real(value))
            # the residual test bounds the multiplicity; `_accept_repeated` decides it (the
            # largest m it accepts)
            multiplicity = kerr_root_multiplicity_at(
                a, energy, lz, q, radius; atol=atol, rtol=rtol)
            accepted = nothing
            while multiplicity >= 2
                accepted = _accept_repeated(a, energy, lz, q, radius, multiplicity,
                    raw_roots, rplus, derivative_order)
                accepted === nothing || break
                multiplicity -= 1
            end
            multiplicity >= max(2, derivative_order + 1) || continue
            radius, reading = accepted
            push!(candidates, (
                radius=radius,
                multiplicity=multiplicity,
                reading=reading,
                source=Symbol("derivative_order_$(derivative_order)"),
                order=derivative_order,
            ))
        end
        degree(derivative_polynomial) > 0 || break
        derivative_polynomial = derivative(derivative_polynomial)
    end

    sort!(candidates; by=item -> item.radius)
    unique_candidates = NamedTuple[]
    for candidate in candidates
        index = findfirst(item -> _root_close(
            item.radius, candidate.radius; atol=10 * atol, rtol=10 * rtol),
            unique_candidates)
        if index === nothing
            push!(unique_candidates, candidate)
        elseif candidate.multiplicity > unique_candidates[index].multiplicity ||
                (candidate.multiplicity == unique_candidates[index].multiplicity &&
                 candidate.order == candidate.multiplicity - 1)
            # a root of multiplicity m is a *simple* root of R^(m-1): located there it is
            # accurate to eps, versus eps^(1/2) as a double root of R^(m-2), etc.
            unique_candidates[index] = candidate
        end
    end
    if !isempty(unique_candidates)
        maximum_multiplicity = maximum(item.multiplicity for item in unique_candidates)
        if maximum_multiplicity >= 3
            best = first(filter(
                item -> item.multiplicity == maximum_multiplicity,
                unique_candidates,
            ))
            return [best], diagnostics
        end
    end
    return unique_candidates, diagnostics
end

"""
    kerr_geo_root_structure(a, E, Lz, Q; atol=1e-12, rtol=1e-12)

Roots of the radial potential R(r): the raw complex roots, the real roots with their
multiplicities and R, R′, R″, R‴ residuals, and their split into roots below, on and
outside the outer horizon r₊. A root has multiplicity m > 1 when R and its first m − 1
derivatives vanish there to the scaled tolerance (`kerr_root_multiplicity_at`); it is then
located as a simple root of R^(m−1). Each root records its `reading`: `:exact`,
`:within_input_ulp` (outside the horizon, constants within one ulp per component of the
repeated-root manifold; inside, a cluster whose replacement changes R(r₊) by less than
the input reach). Nonexact repeated readings describe the repeated-root model, not a
strict solution of the exact supplied constants. All roots are refined together in twice
the working precision (double-double for Float64); distinct roots are not merged on
proximity. The default tolerances are 1e-12 in Float64, carried to other floating-point
types by `_tol`.
"""
function kerr_geo_root_structure(a::Real, energy::Real, lz::Real, q::Real;
        atol::Real=_root_atol(_float_type(a, energy, lz, q)),
        rtol::Real=_root_rtol(_float_type(a, energy, lz, q)))
    horizons = kerr_horizons(a)
    coefficients = kerr_radial_coefficients(a, energy, lz, q)
    polynomial = kerr_radial_polynomial(a, energy, lz, q)
    shifted=_radial_shifted_coefficients(a,energy,lz,q)
    evaluator=_radial_root_evaluator(a,energy,lz,q)
    # starts rotated off the real axis: from real starts the iteration on a real polynomial
    # stays real and cannot reach a complex pair that the companion matrix returned as two
    # close real roots; real roots return to the axis
    # (companion roots of the coefficients rounded to Float64, refined in the working precision)
    T = _float_type(a, energy, lz, q)
    raw_roots = _refine_roots(evaluator, Complex{T}[(1 + Complex{T}(z)) * complex(1, T(2)^-26)
        for z in roots(Polynomial(Float64.(collect(shifted))))])
    repeated, near_repeated = _repeated_root_candidates(
        a, energy, lz, q, polynomial, raw_roots, horizons.rplus; atol=atol, rtol=rtol)

    candidates = NamedTuple[]
    consumed_raw_roots = falses(length(raw_roots))
    # each repeated root takes the m estimates nearest it; one whose nearest estimates were
    # already taken by another repeated root describes the same cluster and is not kept
    for item in sort(repeated; by=item -> -item.multiplicity)
        nearest = sort(eachindex(raw_roots); by=index -> abs(raw_roots[index] - item.radius))
        own = nearest[1:min(item.multiplicity, length(nearest))]
        any(index -> consumed_raw_roots[index], own) && continue
        push!(candidates, item)
        consumed_raw_roots[own] .= true
    end
    # Non-real roots of the real polynomial come in conjugate pairs: when the unconsumed
    # non-real estimates are odd in number, the one nearest the real axis is real.
    leftover = [i for i in eachindex(raw_roots) if !consumed_raw_roots[i]]
    nonreal(i) = abs(imag(raw_roots[i])) > atol + rtol * max(1.0, abs(real(raw_roots[i])))
    unpaired = isodd(count(nonreal, leftover)) ?
        argmin(i -> abs(imag(raw_roots[i])), filter(nonreal, leftover)) : 0
    for (index, value) in pairs(raw_roots)
        consumed_raw_roots[index] && continue
        (!nonreal(index) || index == unpaired) || continue
        radius = _polish_root(coefficients,real(value);evaluator=evaluator)
        push!(candidates, (radius=radius, multiplicity=1, source=:raw_root, reading=:exact))
    end

    # every repeated root consumed its own estimates, so the remaining ones are distinct
    # roots however close they lie; none is merged by proximity
    sort!(candidates; by=item -> item.radius)
    real_roots = NamedTuple[]
    for candidate in candidates
        derivatives = kerr_radial_derivatives(a, energy, lz, q, candidate.radius)
        push!(real_roots, (
            radius=candidate.radius,
            multiplicity=candidate.multiplicity,
            source=candidate.source,
            reading=candidate.reading,
            residuals=(
                R=derivatives.R,
                R1=derivatives.R1,
                R2=derivatives.R2,
                R3=derivatives.R3,
            ),
        ))
    end

    horizon_tolerance = atol + rtol * max(1.0, abs(horizons.rplus))
    near(item) = abs(item.radius - horizons.rplus) <= horizon_tolerance
    # a root within rounding of r₊ is the horizon root when P(r₊) = 0; otherwise R(r₊) = P(r₊)² > 0
    # and that root bounds the allowed sliver just outside the horizon (width P²/(Δ'(r₊)(r₊² + K)))
    on_horizon = _horizon_root(a, energy, lz)
    below = [item for item in real_roots if item.radius < horizons.rplus - horizon_tolerance]
    coincident = on_horizon ? [item for item in real_roots if near(item)] : NamedTuple[]
    exterior = [item for item in real_roots if item.radius > horizons.rplus + horizon_tolerance ||
        (!on_horizon && near(item))]
    real_degree = sum((item.multiplicity for item in real_roots); init=0)

    return (
        polynomial=polynomial,
        degree=degree(polynomial),
        raw_roots=raw_roots,
        real_roots=Tuple(real_roots),
        below_horizon=Tuple(below),
        horizon_coincident=Tuple(coincident),
        exterior=Tuple(exterior),
        complex_root_count=max(0, degree(polynomial) - real_degree),
        multiplicity_sum=real_degree,
        root_tolerance=(atol=atol, rtol=rtol),
        near_repeated=Tuple(near_repeated),
    )
end

# R(r) of a member. When every root is read exactly, R of the given constants evaluated in
# double-double: it keeps its digits next to a turning point, where the factored form loses
# them to the rounding of r − root. A repeated root read from the constants (within one input
# ulp, or an interior cluster merged below rounding) defines the member's model, so then
# R = c ∏(r − xᵢ) over the classified roots (complex ones as (r − ρ)² + η²). So does the
# APEX turning-point geometry (`_apex_root_geometry`), whose roots define the Stable member's
# model. Members' dr/dλ, residuals and radial tables use it.
function _radial_potential_from_roots(a, energy, lz, q, structure)
    # all roots read exactly: R of the given constants, in double-double
    if all(item -> get(item, :reading, :exact) === :exact, structure.real_roots)
        evaluator = _radial_root_evaluator(a, energy, lz, q)
        return r -> first(evaluator(r))
    end
    lead = kerr_radial_coefficients(a, energy, lz, q)[structure.degree + 1]
    reals = map(item -> (item.radius, item.multiplicity), structure.real_roots)
    pairs = [(real(z), imag(z)) for z in _nonreal_roots(structure) if imag(z) > 0]
    return function (r)
        value = lead
        for (x, k) in reals
            value *= (r - x)^k
        end
        for (ρ, η) in pairs
            value *= (r - ρ)^2 + η^2
        end
        return value
    end
end

_potential_scale(a, energy, lz, q, r) =
    max(1.0, _derivative_scales(kerr_radial_coefficients(a, energy, lz, q), r)[1])

function _interval_probe(lower, upper, rplus)
    if isfinite(lower) && isfinite(upper)
        return lower + (upper - lower) / 2
    elseif isfinite(lower)
        return lower + max(1.0, abs(lower), abs(rplus))
    end
    return rplus + max(1.0, abs(rplus))
end

function _allowed_intervals(a, energy, lz, q, root_structure;
        atol=_root_atol(_float_type(a, energy, lz, q)), rtol=_root_rtol(_float_type(a, energy, lz, q)))
    rplus = kerr_horizons(a).rplus
    exterior = root_structure.exterior
    boundaries = [rplus; [item.radius for item in exterior]; Inf]
    potential = _radial_potential_from_roots(a, energy, lz, q, root_structure)
    intervals = NamedTuple[]
    for index in 1:(length(boundaries) - 1)
        lower = boundaries[index]
        upper = boundaries[index + 1]
        probe = _interval_probe(lower, upper, rplus)
        # Root multiplicities fix the sign between roots, even in a narrow
        # forbidden gap whose expanded polynomial is below a residual tolerance.
        value = potential(probe)
        tolerance = atol + rtol * _potential_scale(a, energy, lz, q, probe)
        lower_root = index == 1 ? nothing : exterior[index - 1]
        upper_root = index > length(exterior) ? nothing : exterior[index]
        push!(intervals, (
            lower=lower,
            upper=upper,
            lower_kind=index == 1 ? :outer_horizon : :radial_root,
            upper_kind=isinf(upper) ? :infinity : :radial_root,
            lower_multiplicity=lower_root === nothing ? 0 : lower_root.multiplicity,
            upper_multiplicity=upper_root === nothing ? 0 : upper_root.multiplicity,
            probe=probe,
            potential=value,
            tolerance=tolerance,
            allowed=value >= 0,
        ))
    end
    return Tuple(intervals)
end

_multiplicities(items) = Tuple(item.multiplicity for item in items)

function _case_ids_for_structure(regime, structure)
    below = structure.below_horizon
    exterior = structure.exterior
    below_mult = _multiplicities(below)
    exterior_mult = _multiplicities(exterior)
    complex_count = structure.complex_root_count
    case_ids = Symbol[]

    if regime === :elliptic
        if below_mult == (1,) && exterior_mult == (1, 1, 1) && complex_count == 0
            append!(case_ids, (:A1, :B1))
        elseif below_mult == (1,) && exterior_mult == (1, 2) && complex_count == 0
            append!(case_ids, (:A2, :B2))
        elseif below_mult == (1,) && exterior_mult == (3,) && complex_count == 0
            append!(case_ids, (:K1, :K2))
        elseif below_mult == (1, 1, 1) && exterior_mult == (1,) && complex_count == 0
            push!(case_ids, :B3)
        elseif below_mult == (1,) && exterior_mult == (1,) && complex_count == 2
            push!(case_ids, :B4)
        elseif below_mult == (1,) && exterior_mult == (2, 1) && complex_count == 0
            append!(case_ids, (:K3, :K4, :K5))
        elseif below_mult == (2, 1) && exterior_mult == (1,) && complex_count == 0
            push!(case_ids, :B7)
        elseif below_mult == (1, 2) && exterior_mult == (1,) && complex_count == 0
            push!(case_ids, :B8)
        elseif below_mult == (3,) && exterior_mult == (1,) && complex_count == 0
            push!(case_ids, :B9)
        end
    elseif regime === :parabolic
        if below_mult == (1,) && isempty(exterior_mult) && complex_count == 2
            push!(case_ids, :C1)
        elseif below_mult == (1, 1, 1) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C2)
        elseif below_mult == (1,) && exterior_mult == (2,) && complex_count == 0
            append!(case_ids, (:K6, :K7, :K8))
        elseif below_mult == (1,) && exterior_mult == (1, 1) && complex_count == 0
            append!(case_ids, (:B5, :D1))
        elseif below_mult == (2, 1) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C6)
        elseif below_mult == (1, 2) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C7)
        elseif below_mult == (3,) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C8)
        end
    elseif regime === :hyperbolic
        if isempty(below_mult) && isempty(exterior_mult) && complex_count == 4
            push!(case_ids, :C5)
        elseif below_mult == (1, 1) && isempty(exterior_mult) && complex_count == 2
            push!(case_ids, :C3)
        elseif below_mult == (1, 1, 1, 1) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C4)
        elseif below_mult == (1, 1) && exterior_mult == (2,) && complex_count == 0
            append!(case_ids, (:K9, :K10, :K11))
        elseif below_mult == (1, 1) && exterior_mult == (1, 1) && complex_count == 0
            append!(case_ids, (:B6, :D2))
        elseif below_mult == (1, 2, 1) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C9)
        elseif below_mult == (1, 1, 2) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C10)
        elseif below_mult == (2,) && isempty(exterior_mult) && complex_count == 2
            push!(case_ids, :C11)
        elseif below_mult == (1, 3) && isempty(exterior_mult) && complex_count == 0
            push!(case_ids, :C12)
        end
    end
    return Tuple(case_ids)
end

function _root_endpoint(item; kind=:radial_root, included=true)
    return KerrGeoRadialEndpoint(kind, item.radius, included, item.multiplicity)
end

_horizon_endpoint(rplus) = KerrGeoRadialEndpoint(:outer_horizon, rplus, false, 0)
_infinity_endpoint(::Type{T}) where {T} = KerrGeoRadialEndpoint(:infinity, T(Inf), false, 0)

# A repeated root reached only asymptotically is an open end of the component.
_asymptotic_root(item, kind) = _root_endpoint(item; kind=kind, included=false)

function _case_bounds(case_id, structure, rplus)
    exterior = structure.exterior
    if case_id === :A1
        return _root_endpoint(exterior[2]), _root_endpoint(exterior[3]), :finite_to_finite
    elseif case_id === :A2
        root = _root_endpoint(exterior[2]; kind=:stable_repeated_root)
        return root, root, :constant_radius
    elseif case_id === :K1
        root = _root_endpoint(exterior[1]; kind=:marginal_repeated_root)
        return root, root, :constant_radius
    elseif case_id === :K2
        return _horizon_endpoint(rplus),
            _asymptotic_root(exterior[1], :marginal_repeated_root), :horizon_to_repeated_root
    elseif case_id in (:K3, :K6, :K9)
        root = _root_endpoint(exterior[case_id === :K9 ? end : 1]; kind=:unstable_repeated_root)
        return root, root, :constant_radius
    elseif case_id === :K4
        return _asymptotic_root(exterior[1], :unstable_repeated_root),
            _root_endpoint(exterior[2]), :repeated_root_to_finite
    elseif case_id in (:K5, :K8, :K11)
        return _horizon_endpoint(rplus),
            _asymptotic_root(exterior[1], :unstable_repeated_root), :horizon_to_repeated_root
    elseif case_id in (:K7, :K10)
        return _asymptotic_root(exterior[end], :unstable_repeated_root), _infinity_endpoint(typeof(rplus)),
            :repeated_root_to_infinity
    elseif case_id in (:B1, :B2, :B3, :B4, :B5, :B6, :B7, :B8, :B9)
        return _horizon_endpoint(rplus), _root_endpoint(exterior[1]), :horizon_to_finite
    elseif case_id in (:C1, :C2, :C3, :C4, :C5, :C6, :C7, :C8, :C9, :C10, :C11, :C12)
        return _horizon_endpoint(rplus), _infinity_endpoint(typeof(rplus)), :horizon_to_infinity
    elseif case_id in (:D1, :D2)
        return _root_endpoint(exterior[end]), _infinity_endpoint(typeof(rplus)), :finite_to_infinity
    end
    error("No component bounds are defined for $(case_id).")
end

"""
    kerr_geo_is_critical(component)

The Critical rule: the component's radial motion sits on, or tends asymptotically to, a
repeated root `rc > r₊` of R that is unstable (a double root with R''(rc) > 0) or marginal
(a triple root). A repeated root on the horizon itself (P(r₊) = 0, the horizon and extremal
tiers) does not count. The rule is the same at |a| = 1, where an exterior unstable or
marginal root is Critical as well (K3, K4 and K5 exist at |a| = 1).
"""
kerr_geo_is_critical(c::KerrGeoRadialComponent) =
    _has_critical_end(c.LowerEndpoint, c.UpperEndpoint, c.Metadata.roots, c.Metadata.horizon)

_has_critical_end(lower, upper, roots, rplus) =
    _is_critical_root(lower, roots, rplus) || _is_critical_root(upper, roots, rplus)

function _is_critical_root(endpoint, roots, rplus)
    isfinite(endpoint.Radius) && endpoint.Radius > rplus || return false
    endpoint.Multiplicity == 3 && return true
    endpoint.Multiplicity == 2 || return false
    root = roots[findfirst(r -> r.radius == endpoint.Radius, roots)]
    return root.residuals.R2 > 0
end

function _polar_sector(candidates, requested)
    if requested !== nothing
        requested in candidates || error(
            "Requested polar sector $(requested) is not among $(collect(candidates)).")
        return requested
    end
    length(candidates) == 1 && return first(candidates)
    # Lz = 0 vortical motion is the same geodesic as the axis crossing one.
    candidates == (:axis_crossing, :vortical) && return :axis_crossing
    return :polar_initial_data_required
end

# The same classification with the polar sector `requested` (the choice `kerr_geo_classify`
# would make for it), without solving the roots again: members with a sector rule of their own
# (Q = 0 scatter, C5) relabel the family's classification instead of classifying twice.
function _with_polar_sector(classification::KerrGeoClassification, requested)
    sector = _polar_sector(classification.PolarMetadata.candidates, requested)
    sector === classification.PolarMetadata.selected && return classification
    return KerrGeoClassification(classification.Parameters, classification.ConstantsOfMotion,
        classification.EnergyRegime, classification.MetricLimit, classification.Roots,
        [_with_polar_sector(component, sector) for component in classification.Components],
        classification.CaseIds, classification.ExcludedCaseIds,
        merge(classification.PolarMetadata, (selected=sector,)), classification.Tags,
        classification.SelectedCase, merge(classification.SelectionHint, (polar_sector=requested,)),
        classification.Status)
end

_with_polar_sector(c::KerrGeoRadialComponent, sector) = KerrGeoRadialComponent(c.CaseId,
    c.BroadClass, c.EnergyRegime, c.LowerEndpoint, c.UpperEndpoint, c.Connectivity,
    c.RadialOrientation, c.FormulaFamily, c.PairedCaseIds, c.FamilyMemberCaseIds, sector, c.Tags,
    c.SupportStatus, c.Metadata)

function _base_tags(a, energy, lz, q, structure)
    tags = Symbol[]
    metric = kerr_metric_limit(a)
    metric === :schwarzschild && push!(tags, :schwarzschild_limit)
    metric === :extremal && push!(tags, :extremal_limit)
    metric === :near_extremal && push!(tags, :near_extremal_limit)
    iszero(q) && push!(tags, :zero_carter_constant)
    _zero_lz(a, energy, lz, q) && push!(tags, :zero_axial_angular_momentum)
    any(item -> item.multiplicity >= 2, structure.real_roots) &&
        push!(tags, :repeated_radial_root)
    !isempty(structure.horizon_coincident) && push!(tags, :horizon_coincident_root)
    kerr_energy_regime(energy) === :parabolic && push!(tags, :parabolic_root_at_infinity)
    _on_axis(a, energy, lz, q) && push!(tags, :axis_constants_require_axis_initial_condition)
    return Tuple(unique(tags))
end

function _case_tags(case_id)
    case_id === :A2 && return (:stable_spherical,)
    case_id === :K1 && return (:marginally_stable_spherical, :isso)
    case_id === :K3 && return (:unstable_spherical,)
    case_id === :K4 && return (:homoclinic, :separatrix)
    case_id === :K5 && return (:whirling, :unstable_repeated_root)
    case_id === :C5 && return (:four_complex_radial_roots, :no_radial_turning_point)
    case_id in (:B7, :B8, :B9, :C6, :C7, :C8, :C9, :C10, :C11, :C12) &&
        return (:interior_repeated_root, :degenerate_elementary_radial_formula)
    return ()
end

function _component_for_case(case_id, structure, rplus, energy_regime, polar_sector, base_tags)
    spec = kerr_geo_case(case_id)
    lower, upper, connectivity = _case_bounds(case_id, structure, rplus)
    tags = Tuple(unique((base_tags..., _case_tags(case_id)...)))
    critical = _has_critical_end(lower, upper, structure.real_roots, rplus)
    return KerrGeoRadialComponent(
        case_id,
        critical ? :critical : spec.BroadClass,
        energy_regime,
        lower,
        upper,
        connectivity,
        spec.RadialOrientation,
        spec.FormulaFamily,
        spec.PairedCaseIds,
        spec.FamilyMemberCaseIds,
        polar_sector,
        tags,
        :classified,
        (
            roots=structure.real_roots,
            structure=structure,
            raw_roots=structure.raw_roots,
            root_multiplicities=_multiplicities(structure.real_roots),
            complex_root_count=structure.complex_root_count,
            horizon=rplus,
            allowed_interval=spec.AllowedInterval,
            past_endpoint=spec.PastEndpoint,
            future_endpoint=spec.FutureEndpoint,
        ),
    )
end

function _component_matches_interval(component, interval;
        atol=_root_atol(typeof(component.LowerEndpoint.Radius)),
        rtol=_root_rtol(typeof(component.LowerEndpoint.Radius)))
    root_kinds = (:radial_root, :repeated_root, :stable_repeated_root,
        :marginal_repeated_root, :unstable_repeated_root)
    lower_match = component.LowerEndpoint.Kind == interval.lower_kind ||
        (component.LowerEndpoint.Kind in root_kinds &&
         interval.lower_kind === :radial_root)
    upper_match = component.UpperEndpoint.Kind == interval.upper_kind ||
        (component.UpperEndpoint.Kind in root_kinds &&
         interval.upper_kind === :radial_root)
    lower_radius_match = isinf(component.LowerEndpoint.Radius) == isinf(interval.lower) &&
        (isinf(interval.lower) || _root_close(component.LowerEndpoint.Radius,
            interval.lower; atol=atol, rtol=rtol))
    upper_radius_match = isinf(component.UpperEndpoint.Radius) == isinf(interval.upper) &&
        (isinf(interval.upper) || _root_close(component.UpperEndpoint.Radius,
            interval.upper; atol=atol, rtol=rtol))
    return lower_match && upper_match && lower_radius_match && upper_radius_match
end

function _unassigned_component(interval, regime, polar_sector, base_tags, structure, rplus)
    lower = interval.lower_kind === :outer_horizon ?
        _horizon_endpoint(interval.lower) :
        KerrGeoRadialEndpoint(:radial_root, interval.lower, true,
            interval.lower_multiplicity)
    upper = interval.upper_kind === :infinity ?
        _infinity_endpoint(typeof(rplus)) :
        KerrGeoRadialEndpoint(:radial_root, interval.upper, true,
            interval.upper_multiplicity)
    broad = _has_critical_end(lower, upper, structure.real_roots, rplus) ? :critical :
        interval.lower_kind === :outer_horizon && interval.upper_kind === :infinity ?
        (regime === :elliptic ? :unassigned : :capture) :
        interval.lower_kind === :outer_horizon ? :plunge :
        interval.upper_kind === :infinity ?
            (regime === :elliptic ? :unassigned : :scatter) : :stable
    connectivity = interval.lower_kind === :outer_horizon && interval.upper_kind === :infinity ?
        :horizon_to_infinity :
        interval.lower_kind === :outer_horizon ? :horizon_to_finite :
        interval.upper_kind === :infinity ? :finite_to_infinity : :finite_to_finite
    orientation = broad === :plunge || broad === :capture ? :inward :
        broad === :scatter ? :inbound_turn_outbound :
        broad === :stable ? :libration : :unassigned
    return KerrGeoRadialComponent(
        nothing,
        broad,
        regime,
        lower,
        upper,
        connectivity,
        orientation,
        :UNASSIGNED,
        (),
        (),
        polar_sector,
        Tuple(unique((base_tags..., :unassigned_case))),
        :unassigned,
        (
            roots=structure.real_roots,
            structure=structure,
            raw_roots=structure.raw_roots,
            root_multiplicities=_multiplicities(structure.real_roots),
            complex_root_count=structure.complex_root_count,
            horizon=rplus,
            allowed_interval=(interval.lower, interval.upper),
            past_endpoint=:unassigned,
            future_endpoint=:unassigned,
        ),
    )
end

function _contains_initial_radius(component, radius;
        atol=_root_atol(typeof(component.LowerEndpoint.Radius)),
        rtol=_root_rtol(typeof(component.LowerEndpoint.Radius)))
    lower = component.LowerEndpoint.Radius
    upper = component.UpperEndpoint.Radius
    lower_ok = radius > lower || (component.LowerEndpoint.Included &&
        _root_close(radius, lower; atol=atol, rtol=rtol))
    upper_ok = isinf(upper) || radius < upper || (component.UpperEndpoint.Included &&
        _root_close(radius, upper; atol=atol, rtol=rtol))
    return lower_ok && upper_ok
end

function _radial_sign_matches(component, radial_sign)
    radial_sign === nothing && return true
    radial_sign in (:inward, -1) && return component.RadialOrientation in
        (:inward, :inward_asymptotic, :inbound_turn_outbound, :libration)
    radial_sign in (:outward, 1) && return component.RadialOrientation in
        (:inbound_turn_outbound, :libration, :outward)
    radial_sign in (:zero, 0) && return component.RadialOrientation === :constant_radius
    error("radial_sign must be :inward, :outward, :zero, -1, 0, or 1.")
end

function _endpoint_intent_matches(component, endpoint_intent)
    endpoint_intent === nothing && return true
    endpoint_intent === component.BroadClass && return true
    kinds = (component.LowerEndpoint.Kind, component.UpperEndpoint.Kind)
    endpoint_intent === :horizon && return :outer_horizon in kinds
    endpoint_intent === :infinity && return :infinity in kinds
    endpoint_intent === :finite && return !(:infinity in kinds)
    return false
end

"""
    kerr_geo_select_component(components or classification; case_id, initial_radius,
        radial_sign, endpoint_intent)

The one component that matches every given criterion: `case_id`; `initial_radius` inside its
radial range; `radial_sign` (`:inward`/`-1`, `:outward`/`1`, `:zero`/`0`) allowed by its radial
orientation; `endpoint_intent` (a broad class Symbol, `:horizon`, `:infinity` or `:finite`).
At least one criterion is required; no match or more than one match is an error.
"""
function kerr_geo_select_component(components::AbstractVector;
        case_id=nothing,
        initial_radius=nothing,
        radial_sign=nothing,
        endpoint_intent=nothing)
    any(value -> value !== nothing,
        (case_id, initial_radius, radial_sign, endpoint_intent)) ||
        error("Component selection requires case_id, initial_radius, radial_sign, or endpoint_intent.")
    candidates = collect(components)
    case_id !== nothing && filter!(component -> component.CaseId === case_id, candidates)
    initial_radius !== nothing && filter!(component ->
        _contains_initial_radius(component, float(initial_radius)), candidates)
    radial_sign !== nothing && filter!(component ->
        _radial_sign_matches(component, radial_sign), candidates)
    endpoint_intent !== nothing && filter!(component ->
        _endpoint_intent_matches(component, endpoint_intent), candidates)

    available = [component.CaseId === nothing ? :UNASSIGNED : component.CaseId
        for component in components]
    isempty(candidates) && error(
        "No component matches the requested selection. Available components: $(available).")
    length(candidates) == 1 && return first(candidates)
    matches = [component.CaseId === nothing ? :UNASSIGNED : component.CaseId
        for component in candidates]
    error("Ambiguous component selection matches $(matches). Add case_id or a more specific endpoint intent.")
end

kerr_geo_select_component(classification::KerrGeoClassification; kwargs...) =
    kerr_geo_select_component(classification.Components; kwargs...)

function _classify(a::Real, energy::Real, lz::Real, q::Real;
        polar_sector=nothing,
        case_id=nothing,
        initial_radius=nothing,
        radial_sign=nothing,
        endpoint_intent=nothing,
        atol::Real=_root_atol(_float_type(a, energy, lz, q)),
        rtol::Real=_root_rtol(_float_type(a, energy, lz, q)),
        structure=nothing)
    all(isfinite, (a, energy, lz, q)) || throw(DomainError(
        (a, energy, lz, q), "a, E, Lz and Q must be finite."))
    metric = kerr_metric_limit(a)
    horizons = kerr_horizons(a)
    regime = kerr_energy_regime(energy)
    # the roots of the constants, unless the input supplies its own geometry (APEX turning points)
    structure === nothing && (structure = kerr_geo_root_structure(
        a, energy, lz, q; atol=atol, rtol=rtol))
    intervals = _allowed_intervals(
        a, energy, lz, q, structure; atol=atol, rtol=rtol)
    polar = kerr_polar_admissibility(a, energy, lz, q)
    polar_candidates = kerr_polar_sector_candidates(a, energy, lz, q)
    selected_polar = _polar_sector(polar_candidates, polar_sector)
    horizon_momentum = kerr_radial_momentum(a, energy, lz, horizons.rplus)
    future_horizon = _horizon_root(a, energy, lz) ? :horizon_root :
        horizon_momentum > 0 ? :future_directed : :past_directed
    candidate_ids = _case_ids_for_structure(regime, structure)
    base_tags = _base_tags(a, energy, lz, q, structure)

    components = KerrGeoRadialComponent[]
    # E < 0 (Class N) and horizon-root constants have their own constructors (`kerr_geodesic`)
    built_elsewhere = energy < 0 || !isempty(structure.horizon_coincident)
    if polar.admissible && !built_elsewhere
        for id in candidate_ids
            component = _component_for_case(
                id, structure, horizons.rplus, regime, selected_polar, base_tags)
            if component.LowerEndpoint.Kind === :outer_horizon &&
                    future_horizon !== :future_directed
                continue
            end
            push!(components, component)
        end

        for interval in intervals
            interval.allowed || continue
            any(component -> _component_matches_interval(
                component, interval; atol=atol, rtol=rtol), components) && continue
            if interval.lower_kind === :outer_horizon &&
                    future_horizon !== :future_directed
                continue
            end
            push!(components, _unassigned_component(
                interval, regime, selected_polar, base_tags, structure, horizons.rplus))
        end
    end

    selection_hint = (
        case_id=case_id,
        initial_radius=initial_radius,
        radial_sign=radial_sign,
        endpoint_intent=endpoint_intent,
        polar_sector=polar_sector,
    )
    has_selection = any(value -> value !== nothing,
        (case_id, initial_radius, radial_sign, endpoint_intent))
    selected = has_selection ? kerr_geo_select_component(
        components;
        case_id=case_id,
        initial_radius=initial_radius,
        radial_sign=radial_sign,
        endpoint_intent=endpoint_intent,
    ) : nothing
    case_ids = Tuple(component.CaseId for component in components
        if component.CaseId !== nothing)
    unassigned_count = count(component -> component.CaseId === nothing, components)
    reason = !polar.admissible ? :polar_motion_inadmissible :
        energy < 0 ? :trapped_class_n :
        !isempty(structure.horizon_coincident) ? :horizon_coincident_root :
        isempty(components) ? :no_future_directed_exterior_component :
        unassigned_count > 0 ? :classified_with_unassigned_component : :classified

    return KerrGeoClassification(
        (a=float(a), E=float(energy), Lz=float(lz), Q=float(q)),
        (E=float(energy), Lz=float(lz), Q=float(q)),
        regime,
        metric,
        structure.raw_roots,
        components,
        case_ids,
        (),
        (
            admissibility=polar,
            candidates=polar_candidates,
            selected=selected_polar,
        ),
        base_tags,
        selected === nothing ? nothing : selected.CaseId,
        selection_hint,
        (
            classified=reason in (:classified, :classified_with_unassigned_component),
            reason=reason,
            energy_sign=kerr_energy_sign(energy),
            future_horizon_condition=future_horizon,
            horizon_momentum=horizon_momentum,
            root_structure=structure,
            repeated_root_reading=any(x -> x.multiplicity > 1 && x.reading !== :exact,
                structure.real_roots) ? :within_input_ulp : :exact,
            allowed_intervals=intervals,
            unassigned_component_count=unassigned_count,
            selected_component=selected,
        ),
    )
end

"""
    kerr_geo_classify(a, E, Lz, Q; polar_sector, case_id, initial_radius, radial_sign,
        endpoint_intent, atol, rtol) -> KerrGeoClassification

Roots of R(r), the allowed exterior radial intervals and the `KerrGeoRadialComponent`s, each
with case ID, broad class and polar sector: one per allowed interval (an interval ending on the
horizon only when P(r₊) > 0) and one per constant-radius orbit on a repeated root. The
selection keywords (as in `kerr_geo_select_component`) set `SelectedCase` and keep the other
components. Constants with E < 0 (Class N), with a root of R on the outer horizon, or with
inadmissible polar motion have no components here, and `Status.reason` says which;
`kerr_geodesic` builds the Class N and horizon-root members.
"""
function kerr_geo_classify(a::Real, energy::Real, lz::Real, q::Real; kwargs...)
    return _classify(a, energy, lz, q; kwargs...)
end

"""
    kerr_geo_components(a, E, Lz, Q; kwargs...)

`kerr_geo_classify(a, E, Lz, Q; kwargs...).Components`: the classified components; selection
keywords do not remove any.
"""
function kerr_geo_components(a::Real, energy::Real, lz::Real, q::Real; kwargs...)
    return kerr_geo_classify(a, energy, lz, q; kwargs...).Components
end
