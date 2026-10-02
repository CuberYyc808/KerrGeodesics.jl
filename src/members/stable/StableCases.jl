# Class A (Stable) helpers: the APEX inclination x from the constants, the constants residual,
# selection of the Stable component and its stability metadata.

# (E < 1 and Q ≥ 0: the discriminant (Q − c)² + 2(Q + c)Lz² + Lz⁴ is never negative)
function _polar_turning_cosine_squared(a, energy, lz, q)
    _zero_lz(a, energy, lz, q) && return 1.0
    roots = _polar_quadratic_roots(a, energy, lz, q)
    iszero(roots.c) && iszero(q + lz^2) && error(
        "The polar turning point is not isolated: a²(E² − 1) and Q + Lz² both vanish.")
    admissible = [clamp(value, 0.0, 1.0) for value in (roots.u_small, roots.u_big)
        if -POLAR_ROOT_SLACK <= value <= 1 + POLAR_ROOT_SLACK]
    isempty(admissible) && error(
        "Polar APEX inversion found no turning point in cos(theta)^2 in [0,1].")
    return maximum(admissible)
end

function _apex_x(a, energy, lz, q)
    _zero_lz(a, energy, lz, q) && return 0.0
    uturn = _polar_turning_cosine_squared(a, energy, lz, q)
    # x² = 1 − z²_turn, formed without cancellation (near-polar orbits: x → 0)
    return sign(lz) * sqrt(max(0.0, _polar_one_minus_root(a, energy, lz, q, uturn)))
end

function _constants_residual(input, reconstructed)
    residuals = (
        E=reconstructed.E - input.E,
        Lz=reconstructed.Lz - input.Lz,
        Q=reconstructed.Q - input.Q,
    )
    scaled = (
        E=abs(residuals.E) / max(1.0, abs(input.E)),
        Lz=abs(residuals.Lz) / max(1.0, abs(input.Lz)),
        Q=abs(residuals.Q) / max(1.0, abs(input.Q)),
    )
    return (
        absolute=residuals,
        scaled=scaled,
        maximum_scaled=max(scaled.E, scaled.Lz, scaled.Q),
    )
end

function _stable_radial_component(classification, case_id)
    candidates = [component for component in classification.Components
        if component.BroadClass === :stable]
    case_id !== nothing && filter!(component -> component.CaseId === case_id, candidates)
    isempty(candidates) && error(
        "No Class A component is available. Classified cases: $(classification.CaseIds).")
    length(candidates) == 1 || error(
        "More than one Class A component is available; provide case_id. Cases: $([item.CaseId for item in candidates]).")
    return only(candidates)
end


"""
    kerr_geo_stability_metadata(a, E, Lz, Q, component; atol=1e-8)

Shape and radial stability of a Stable component (A1 or A2): `shape` (`:eccentric`,
`:circular` or `:spherical`), `stability` (`:stable` on the double root of A2, where
R''(r) < 0 is checked; not applicable to the eccentric A1), the `radial_derivatives` at that
root, and the `formula_family`.
"""
function kerr_geo_stability_metadata(a::Real, energy::Real, lz::Real, q::Real,
        component::KerrGeoRadialComponent; atol::Real=1.0e-8)
    component.CaseId in (:A1, :A2) || error(
        "Stability metadata requires a Class A component.")
    shape = component.CaseId === :A1 ? :eccentric :
        iszero(q) ? :circular : :spherical
    if component.CaseId === :A1
        stability = :not_applicable_nonconstant_radius
        radial_derivatives = nothing
        limit = :none
        stability_check_passed = true
    else
        radius = component.LowerEndpoint.Radius
        radial_derivatives = kerr_radial_derivatives(a, energy, lz, q, radius)
        stability = :stable
        limit = :none
        scale = max(1.0, abs(radial_derivatives.R3) * max(1.0, abs(radius)))
        stability_check_passed = radial_derivatives.R2 < -atol * scale
        stability_check_passed || error(
            "The repeated root of $(component.CaseId) is not radially stable " *
            "(R'' must be negative).")
    end
    return (
        broad_class=:stable,
        case_id=component.CaseId,
        stability=stability,
        stability_scope=component.CaseId === :A1 ?
            :constant_radius_test_not_applicable : :radial_repeated_root,
        shape=shape,
        limit=limit,
        radial_derivatives=radial_derivatives,
        stability_check_passed=stability_check_passed,
        formula_family=component.FormulaFamily,
    )
end
