# Class A (Stable) helpers: the APEX inclination x from the constants, the constants residual,
# selection of the Stable component and its stability metadata.

function _polar_turning_cosine_squared(a, energy, lz, q; atol=1.0e-11)
    abs(lz) <= atol && return 1.0
    beta = a^2 * _e2m1(energy)
    quadratic = -beta
    linear = beta - q - lz^2
    roots_u = Float64[]
    if abs(quadratic) <= atol
        abs(linear) <= atol && error(
            "The polar turning point is not isolated: a²(E² − 1) and Q + Lz² both vanish.")
        push!(roots_u, -q / linear)
    else
        discriminant = linear^2 - 4 * quadratic * q
        discriminant >= -atol || error(
            "Polar APEX inversion has no real turning point.")
        root = sqrt(max(0.0, discriminant))
        push!(roots_u, (-linear - root) / (2 * quadratic))
        push!(roots_u, (-linear + root) / (2 * quadratic))
    end
    admissible = [clamp(value, 0.0, 1.0) for value in roots_u
        if -atol <= value <= 1 + atol]
    isempty(admissible) && error(
        "Polar APEX inversion found no turning point in cos(theta)^2 in [0,1].")
    return maximum(admissible)
end

function _apex_x(a, energy, lz, q; atol=1.0e-11)
    abs(lz) <= atol && return 0.0
    uturn = _polar_turning_cosine_squared(a, energy, lz, q; atol=atol)
    # x² = 1 − z²_turn, formed without cancellation (near-polar orbits: x → 0)
    return sign(lz) * sqrt(max(0.0, _polar_one_minus_root(a, energy, lz, q, uturn)))
end

function _constants_residual(input, reconstructed)
    residuals = (
        E=reconstructed["E"] - input.E,
        Lz=reconstructed["Lz"] - input.Lz,
        Q=reconstructed["Q"] - input.Q,
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


function kerr_geo_stability_metadata(a::Real, energy::Real, lz::Real, q::Real,
        component::KerrGeoRadialComponent; atol::Real=1.0e-8)
    component.CaseId in (:A1, :A2) || error(
        "Stability metadata requires a Class A component.")
    equatorial = abs(q) <= atol
    shape = component.CaseId === :A1 ? :eccentric :
        equatorial ? :circular : :spherical
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
