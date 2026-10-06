# Plunge reference API: Jacobi r(λ), θ(λ) of an E < 1 plunge; t, φ from the radial engine
# (one radial period, continued through both horizons) and the polar engine.

function real2_radial_position(absλ, E, roots, Kr, Jr)
    r4, r3, r2, r1 = roots
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    sn2 = _sn(Kr - ξr * absλ, Jr)^2
    return (r3 * (r1 - r2) * sn2 - r2 * (r1 - r3)) /
           ((r1 - r2) * sn2 - (r1 - r3))
end

"""
    generic_plunge_orbit(a, E, L, Q; initPhases=(0.0, 0.0, 0.0, 0.0))

Return the Boyer-Lindquist `t(λ)`, `r(λ)`, `theta(λ)` and `phi(λ)` of an E < 1 Kerr plunge
in the root class of `classify_orbit` (Real1, Real2 or Complex).

r(λ) and θ(λ) are the Jacobi closed forms; t and φ come from the radial spectral engine
(one radial period, outer turning point → inner turning point → outer turning point,
continued through r₊ and r₋ as principal values) and the polar engine.
"""
function generic_plunge_orbit(a, E, L, Q; initPhases = (0.0, 0.0, 0.0, 0.0), real2_horizon_offset=1e-4)
    orbit = _plunge_orbit(a, E, L, Q; initPhases=initPhases)
    return orbit === nothing ? nothing : [orbit.t, orbit.r, orbit.theta, orbit.phi]
end

# t, r, θ, φ of the plunge and v = t + r* from the radial engine's ingoing regular chart, finite
# and free of the cancellation of t + r* up to the future horizon
function _plunge_orbit(a, E, L, Q; initPhases = (0.0, 0.0, 0.0, 0.0))
    roots, cf = classify_orbit(a, E, L, Q)
    λr0 = initPhases[2]
    if cf == "Real1" || cf == "Real2"
        r4, r3, r2, r1 = roots
        ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
        kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
        Kr, Jr = _K(kr), _jacobi_parameter(kr)
        half = Kr / ξr                             # outer turning point → inner one
        cf == "Real2" && return _generic_plunge_orbit(a, E, L, Q, initPhases, half,
            λ -> real2_radial_position(λ + λr0, E, roots, Kr, Jr))
        return _generic_plunge_orbit(a, E, L, Q, initPhases, half, function (λ)
            sn2 = _sn(ξr * (λ + λr0), Jr)^2
            return (r3 * (r2 - r4) - r2 * (r3 - r4) * sn2) / ((r2 - r4) - (r3 - r4) * sn2)
        end)
    elseif cf == "Complex"
        r1, r2, A, B = roots
        ξr = sqrt((1 - E^2) * A * B)
        kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)
        Jr = _jacobi_parameter(kr)
        return _generic_plunge_orbit(a, E, L, Q, initPhases, 2 * _K(kr) / ξr,
            function (λ)
                sn, cn, _ = _jacobi_sncndn(ξr * (λ + λr0), Jr)
                return (2 * A * B * (r1 + r2) + (A - B) * (A * r2 - B * r1) * sn^2 +
                    2 * A * B * (r1 - r2) * cn) / (4 * A * B + (A - B)^2 * sn^2)
            end)
    end
    @info("generic_plunge_orbit: root class $cf has no E < 1 plunge region.")
    return nothing
end

function _generic_plunge_orbit(a, E, L, Q, initPhases, half, radius)
    λt0, λr0, λθ0, λϕ0 = initPhases
    zm, ξθ, kθ = _plunge_polar_parameters(a, E, L, Q)
    Jθ = _jacobi_parameter(kθ)
    θ(λ) = acos(sqrt(zm) * _sn(ξθ * (λ + λθ0), Jθ))
    # z = √z₋ sn(ξθ(λ + λθ0)): the polar engine's phase is counted from the northern turning point
    polar = _polar_solution(a, E, L, Q, iszero(Q) ? :equatorial : :pendular,
        ξθ * λθ0 - _K(kθ))
    # one radial period containing λ = 0, starting at an outer turning point (λ = −λr0 mod 2·half)
    period = 2half
    start = -λr0 - period * floor(-λr0 / period)
    start > 0 && (start -= period)
    coords = _engine_coordinates(a, E, L, Q, radius, _polar_primitive(polar);
        potential=_coefficient_potential(a, E, L, Q),
        domain=(start, start + period), ends=(:turning, :turning), turn=start + half,
        σ=1.0, period=period, λ_bl=0.0)
    t(λ) = _coords_t(coords, λ) + λt0
    ϕ(λ) = _coords_phi(coords, λ) + λϕ0
    # v − v(start) from the regular chart (t + r*, zero at the outer turning point `start`)
    chart = _regular_chart(coords, -one(start), start)
    v(λ) = chart(λ)[1] + (_coords_t(coords, start) + _rstar_all(a, radius(start))) + λt0
    return (t=t, r=radius, theta=θ, phi=ϕ, v=v)
end
