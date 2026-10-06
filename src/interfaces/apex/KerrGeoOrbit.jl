# APEX reference API: `kerr_geo_orbit`, the closed-form trajectory in (a, p, e, x) of stable
# orbits and of the constant-radius Critical orbits with E < 1 (ISCO/ISSO, unstable circular and
# spherical orbits).

"""
Orbit-type metadata of (a, p, e, x) when `kerr_geo_orbit` builds a Dict for it, otherwise an
ArgumentError naming the orbit's class. Constant-radius Critical orbits with E < 1 are built
like stable ones; their radial phase does not advance, so ϒr = 0 (`_constant_radius_dict!`).
"""
function _apex_orbit_gate(a, p, e, x)
    meta = kerr_geo_orbit_type_metadata(a, p, e, x)
    meta.family == "Stable" && return meta
    if meta.family == "Critical" && meta.shape == "Circular"
        meta.energy_regime == "Elliptic" && return meta
        throw(ArgumentError("kerr_geo_orbit: (a, p, e, x) = $((a, p, e, x)) is an unstable " *
            "constant-radius orbit with E ≥ 1 ($(meta.energy_regime)), outside the range of " *
            "the APEX frequency formulas. Build it with kerr_geodesic(a, p, e, x): Critical " *
            "member K6 (E = 1) or K9 (E > 1)."))
    end
    throw(ArgumentError("kerr_geo_orbit: (a, p, e, x) = $((a, p, e, x)) is a $(meta.family) " *
        "orbit ($(join(meta.labels, ", "))). kerr_geo_orbit builds Stable orbits and " *
        "constant-radius Critical orbits with E < 1; use " *
        "kerr_geodesic(a, p, e, x) for the members of this one."))
end

# A constant-radius Critical orbit: ϒr = 0 (the radial phase does not advance), with its
# stability, and for an unstable one the Mino-time growth rate √(R''(p)/2) of a radial
# perturbation (δr'' = R''(p) δr/2).
function _constant_radius_dict!(orbit, meta)
    orbit["Stability"] = meta.stability
    meta.family == "Critical" || return orbit
    orbit["RadialFrequency"] = 0.0
    orbit["Frequencies"]["ϒr"] = 0.0
    if meta.stability == "Unstable"
        R2 = kerr_radial_derivatives(orbit["a"], orbit["Energy"], orbit["AngularMomentum"],
            orbit["CarterConstant"], meta.effective_p).R2
        orbit["RadialLyapunovExponent"] = sqrt(max(R2, 0.0) / 2)
    end
    return orbit
end

# ---- Kerr circular orbit (Mino time) ----

function kerr_geo_orbit_circular(a::Real, p::Real, e::Real=0.0, x::Real=1.0; initPhases=(0.0, 0.0, 0.0, 0.0))
    # Orbit type
    orbit_type = _apex_orbit_gate(a, p, e, x)
    type = orbit_type.labels

    # Orbital frequencies (Mino time)
    Frequencies = kerr_geo_frequencies(a, p, e, x; Time="Mino")
    ϒt = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒϕ = Frequencies["ϒϕ"] 
    # Constants of motion
    consts = kerr_geo_constants_of_motion(a, p, e, x)
    En, Lz, Q = consts["E"], consts["Lz"], consts["Q"]
    # Radial roots
    r1,r2,r3,r4 = kerr_geo_radial_roots(a, p, e, x; En = En, Lz = Lz, Q = Q)

    # Trajectory functions (broadcastable)
    # Uniform motion at the Mino frequencies (prograde and retrograde, x = ±1).
    t(λ) = initPhases[1] + ϒt * λ
    r(λ) = p  # constant radius
    θ(λ) = oftype(ϒt, π)/2 # equatorial
    ϕ(λ) = initPhases[4] + ϒϕ * λ

    # Four-velocity
    velocity = kerr_geo_four_velocity(a, p, e, x; initPhases=(initPhases[2], initPhases[3]), Covariant=false, Parametrization="Mino")

    # Associate dictionary
    return _constant_radius_dict!(Dict{String,Any}(
        "Parametrization" => "Mino",
        "Energy" => En,
        "AngularMomentum" => Lz,
        "CarterConstant" => Q,
        "ConstantsOfMotion" => consts,
        "RadialRoots" => [r1, r2, r3, r4],
        "RadialFrequency" => ϒr,
        "PolarFrequency" => ϒθ,
        "AzimuthalFrequency" => ϒϕ,
        "Frequencies" => Dict("ϒt" => ϒt, "ϒr" => ϒr, "ϒθ" => ϒθ, "ϒϕ" => ϒϕ),
        "Trajectory" => [t, r, θ, ϕ],
        "FourVelocity" => velocity,
        "CrossFunction" => nothing,
        "DerivativesCrossFunction" => nothing,
        "a" => a,
        "p" => p,
        "e" => e,
        "Cosθ_inc" => x,
        "Type" => type,
        "InitialPhases" => initPhases
    ), orbit_type)
end


function kerr_geo_orbit_generic(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    # Orbit type
    orbit_type = _apex_orbit_gate(a, p, e, x)
    type = orbit_type.labels
    # Get constants of motion: Energy, angular momentum, Carter constant
    consts = kerr_geo_constants_of_motion(a, p, e, x)
    En, Lz, Q = consts["E"], consts["Lz"], consts["Q"]

    # Get Mino-time fundamental frequencies
    Frequencies = kerr_geo_frequencies(a, p, e, x; Time="Mino")
    ϒt = Frequencies["ϒt"]
    ϒr = Frequencies["ϒr"]
    ϒθ = Frequencies["ϒθ"]
    ϒϕ = Frequencies["ϒϕ"]

    # Radial and polar roots
    r1,r2,r3,r4 = kerr_geo_radial_roots(a, p, e, x; En, Lz, Q)
    zp, zm = kerr_geo_polar_roots(a, p, e, x)

    # Jacobi elliptic modulus for radial and polar motion
    kr = ((r1-r2)/(r1-r3)) * ((r3-r4)/(r2-r4))
    kθ = a^2 * (1-En^2) * (zm/zp)^2

    # Horizon radii
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)

    # Elliptic Pi parameters for radial motion
    hr = (r1-r2)/(r1-r3)
    hp = ((r1-r2)*(r3-rp))/((r1-r3)*(r2-rp))
    hm = ((r1-r2)*(r3-rm))/((r1-r3)*(r2-rm))

    # complete integrals and Jacobi parameter records, formed once for the orbit
    Kr, Kθ = _K(kr), _K(kθ)
    Jr, Jθ = _jacobi_parameter(kr), _jacobi_parameter(kθ)
    halfπ = oftype(Kθ, π) / 2

    # Radial JacobiSN mapping
    rq(qr) = (r3*(r1 - r2) * _sn(Kr/π * qr, Jr)^2 - r2*(r1-r3)) /
            ((r1-r2) * _sn(Kr/π * qr, Jr)^2 - (r1-r3))

    # Polar JacobiSN mapping
    zq(qθ) = zm * _sn(Kθ * 2/π * (qθ + halfπ), Jθ)

    # Radial and polar Jacobi amplitudes
    ψr(qr) = _am(Kr/π * qr, Jr)
    ψθ(qθ) = _am(Kθ*2/π*(qθ+halfπ), Jθ)
    dψr(qr) = _dn(Kr/π * qr, Jr) * Kr / π
    dψθ(qθ) = 2 * _dn(Kθ*2/π*(qθ+halfπ), Jθ) * Kθ / π

    # t and phi increments due to radial motion

    function elliptic_pi_r(h, k, qr)
        complete = _Pi(h, k)
        ψ = ψr(qr)
        period = fld(ψ, oftype(ψ, π))   # floor: remainder in [0, π) also for λ < 0
        remainder = ψ - period * oftype(ψ, π)
        if remainder <= oftype(ψ, π) / 2
            incomplete = _Pi(h, remainder, k)
        else
            remainder = oftype(ψ, π) - remainder
            incomplete = 2 * complete - _Pi(h, remainder, k)
        end
        incomplete += period * complete * 2
        return complete * qr / π - incomplete
    end

    function elliptic_pi_θ(h, k, qθ)
        complete = _Pi(h, k)
        ψ = ψθ(qθ)
        period = fld(ψ, oftype(ψ, π))   # floor: remainder in [0, π) also for λ < 0
        remainder = ψ - period * oftype(ψ, π)
        if remainder <= oftype(ψ, π) / 2
            incomplete = _Pi(h, remainder, k)
        else
            remainder = oftype(ψ, π) - remainder
            incomplete = 2 * complete - _Pi(h, remainder, k)
        end
        incomplete += period * complete * 2
        return complete*2*((qθ+halfπ)/π) - incomplete
    end

    # Horizon terms enter as [F(r+) - F(r-)]/(r+ - r-). At |a| -> 1 both horizons merge and
    # this becomes F'(1); _horizon_divdiff takes that limit (Richardson-extrapolated).
    hρ(ρ) = ((r1 - r2) * (r3 - ρ)) / ((r1 - r3) * (r2 - ρ))
    Fsq(ρ) = (r2 - ρ) * (r3 - ρ)
    sqrtfac = sqrt(max(0.0, (1 - En^2) * (r1 - r3) * (r2 - r4)))
    # Spherical orbits (e = 0): r is constant, so the radial oscillation terms vanish
    # (the elliptic forms degenerate at a repeated or triple root).
    spherical = isapprox(e, 0; atol=_tol(typeof(En), 1e-12))
    o = zero(En)

    function Δtr(qr)
        spherical && return o
        prefac = -En / sqrtfac
        term1 = 4 * (r2 - r3) * elliptic_pi_r(hr, kr, qr)
        F(ρ) = (-2a^2 + ρ * (4 - a * Lz / En)) / Fsq(ρ) * elliptic_pi_r(hρ(ρ), kr, qr)
        term2 = -4 * (r2 - r3) * _horizon_divdiff(F, rp, rm)
        term3 = (r2 - r3) * (r1 + r2 + r3 + r4) * elliptic_pi_r(hr, kr, qr)
        term4 = (r1 - r3) * (r2 - r4) * (_E(kr) * qr / π - _E(ψr(qr), kr) +
            hr * (sin(ψr(qr)) * cos(ψr(qr)) * sqrt(1 - kr * sin(ψr(qr))^2)) / (1 - hr * sin(ψr(qr))^2))
        return prefac * (term1 + term2 + term3 + term4)
    end

    function Δϕr(qr)
        spherical && return o
        G(ρ) = (2ρ - a * Lz / En) / Fsq(ρ) * elliptic_pi_r(hρ(ρ), kr, qr)
        return 2a * En * (r2 - r3) / sqrtfac * _horizon_divdiff(G, rp, rm)
    end

    # t and phi increments due to polar motion
    # En zp/(1-En^2) (E_c x - E(ψ)) = En a^2 zm^2/zp (D(ψ) - D_c x), x = 2(qθ+π/2)/π, since
    # F(ψ) = K x; the D form has no 1/(1 - En^2) (wide orbits, En -> 1).
    Dθc = _D(kθ)
    Δtθ(qθ) = En * a^2 * zm^2 / zp * (_D(ψθ(qθ), kθ) - Dθc * 2 * ((qθ+halfπ)/π))
    Δϕθ(qθ) = iszero(Lz) ? o : -Lz/zp * elliptic_pi_θ(zm^2, kθ, qθ)   # Lz = 0: polar orbit, zm = 1

    ellip_diff_pi_r(h, k, q) = _Pi(h, k)/π -
        dψr(q)/((1 - h * sin(ψr(q))^2) * sqrt(1 - k * sin(ψr(q))^2))
    function dtr(qr)
        spherical && return o
        ψ = ψr(qr); sn2 = sin(ψ)^2; cn2 = cos(ψ)^2; w = sqrt(1 - kr * sn2); dψ = dψr(qr)
        regular = 4 * (r2 - r3) * ellip_diff_pi_r(hr, kr, qr) +
            (r2 - r3) * (r1 + r2 + r3 + r4) * ellip_diff_pi_r(hr, kr, qr) +
            (r1 - r3) * (r2 - r4) * (_E(kr)/π - (hr * kr * cn2 * sn2 * dψ)/((1 - hr * sn2) * w) -
                w * dψ + (2 * hr^2 * cn2 * sn2 * w * dψ)/(1 - hr * sn2)^2 +
                (hr * cn2 * w * dψ)/(1 - hr * sn2) - (hr * sn2 * w * dψ)/(1 - hr * sn2))
        F(ρ) = (-2a^2 + (4 - a * Lz / En) * ρ) / Fsq(ρ) * ellip_diff_pi_r(hρ(ρ), kr, qr)
        return -En * (regular - 4 * (r2 - r3) * _horizon_divdiff(F, rp, rm)) / sqrtfac
    end

    function dϕr(qr)
        spherical && return o
        G(ρ) = (2ρ - a * Lz / En) / Fsq(ρ) * ellip_diff_pi_r(hρ(ρ), kr, qr)
        return 2a * En * (r2 - r3) / sqrtfac * _horizon_divdiff(G, rp, rm)
    end

    dtθ(qθ) = En * a^2 * zm^2 / zp * (sin(ψθ(qθ))^2 * dψθ(qθ) / sqrt(1 - kθ * sin(ψθ(qθ))^2) -
        2 * Dθc / π)
    dϕθ(qθ) = iszero(Lz) ? o : (-(2 * Lz * _Pi(zm^2, kθ) / π) + (Lz * dψθ(qθ)) / (sqrt(1 - kθ * sin(ψθ(qθ))^2) * (1 - zm^2 * sin(ψθ(qθ))^2))) / zp

    qt0, qr0, qθ0, qϕ0 = initPhases

    # Total trajectory functions
    t(λ) = qt0 + ϒt * λ + Δtr(qr0 + ϒr * λ) + Δtθ(qθ0 + ϒθ * λ) - Δtr(qr0) - Δtθ(qθ0)
    r(λ) = rq(qr0 + ϒr * λ)
    θ(λ) = acos(zq(qθ0 + ϒθ * λ))
    ϕ(λ) = qϕ0 + ϒϕ * λ + Δϕr(qr0 + ϒr * λ) + Δϕθ(qθ0 + ϒθ * λ) - Δϕr(qr0) - Δϕθ(qθ0)

    # Four-velocity
    velocity = kerr_geo_four_velocity(a, p, e, x; initPhases=(initPhases[2], initPhases[3]), Covariant=false, Parametrization="Mino")

    return _constant_radius_dict!(Dict{String,Any}(
        "a" => a,
        "p" => p,
        "e" => e,
        "Cosθ_inc" => x,
        "Parametrization" => "Mino",
        "Energy" => En,
        "AngularMomentum" => Lz,
        "CarterConstant" => Q,
        "ConstantsOfMotion" => consts,
        "RadialRoots" => [r1,r2,r3,r4],
        "RadialFrequency" => ϒr,
        "PolarFrequency" => ϒθ,
        "AzimuthalFrequency" => ϒϕ,
        "Frequencies" => Dict("ϒt" => ϒt, "ϒr" => ϒr, "ϒθ" => ϒθ, "ϒϕ" => ϒϕ),
        "Trajectory" => [t,r,θ,ϕ],
        "CrossFunction" => [Δtr, Δtθ, Δϕr, Δϕθ],
        "DerivativesCrossFunction" => [dtr, dtθ, dϕr, dϕθ],
        "FourVelocity" => velocity,
        "Type" => type,
        "InitialPhases" => initPhases
    ), orbit_type)
end

"""
    kerr_geo_orbit(a, p, e, x; initPhases=(0.0, 0.0, 0.0, 0.0), precision=nothing)

The stable orbit with APEX parameters `(a, p, e, x)`, or the constant-radius Critical orbit
with E < 1 (the ISCO or ISSO, an unstable circular or spherical orbit), in closed form as a
`Dict`:

- `"Trajectory"`: `[t, r, θ, ϕ]`, functions of Mino time λ; `"FourVelocity"`: `[uᵗ, uʳ, uᶿ, uᵠ]`;
- `"Frequencies"`: the Mino frequencies `"ϒt"`, `"ϒr"`, `"ϒθ"`, `"ϒϕ"`, also stored as
  `"RadialFrequency"`, `"PolarFrequency"`, `"AzimuthalFrequency"`;
- `"CrossFunction"` and `"DerivativesCrossFunction"`: the oscillating parts `[Δtr, Δtθ, Δϕr,
  Δϕθ]` of t and ϕ and their derivatives (`nothing` on a circular equatorial orbit);
- `"Energy"`, `"AngularMomentum"`, `"CarterConstant"`, `"ConstantsOfMotion"`, `"RadialRoots"`,
  `"a"`, `"p"`, `"e"`, `"Cosθ_inc"`, `"Type"`, `"InitialPhases"`, `"Parametrization"`;
- `"Stability"`; a constant-radius Critical orbit has `"ϒr" = 0` and, when unstable,
  `"RadialLyapunovExponent"` = √(R''/2), the Mino-time growth rate of a radial perturbation.

`initPhases = (qt0, qr0, qθ0, qϕ0)`. The orbit is computed in the floating-point type of
`(a, p, e, x)`; `precision = p` converts them to `BigFloat` of `p` bits, and the returned
functions evaluate at that precision. Any other orbit throws an `ArgumentError`;
`kerr_geodesic(a, p, e, x)` builds every class, including unbound orbits (e > 1).
"""
function kerr_geo_orbit(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0),
        precision=nothing)
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geo_orbit(BigFloat(a), BigFloat(p), BigFloat(e), BigFloat(x); initPhases=initPhases)
    end
    T = _float_type(a, p, e, x)
    a, p, e, x = T(a), T(p), T(e), T(x)
    e > 1 && throw(ArgumentError("kerr_geo_orbit: e = $e > 1 is an unbound orbit; the bound-orbit " *
        "formulas hold for e ≤ 1. Use kerr_geodesic(a, p, e, x).Scatter."))
    prec = _input_precision(a, p, e, x)
    orbit = _with_precision(T, prec) do
        if isapprox(e, 0; atol = _tol(T, 1e-12)) && isapprox(abs(x), 1; atol = _tol(T, 1e-12))
            kerr_geo_orbit_circular(a, p, e, x; initPhases = initPhases)
        else
            kerr_geo_orbit_generic(a, p, e, x; initPhases = initPhases)
        end
    end
    # functions of λ evaluate at the precision the orbit was built with
    for key in ("Trajectory", "FourVelocity", "CrossFunction", "DerivativesCrossFunction")
        orbit[key] = _precision_wrap(T, prec, orbit[key])
    end
    return orbit
end
