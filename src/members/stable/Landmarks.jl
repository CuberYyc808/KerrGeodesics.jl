# Spherical-orbit landmarks: ISCO, ISSO, photon sphere, IBSO, separatrix, and the orbit-type
# metadata of an APEX (p, e, x) orbit.

# Innermost stable circular orbit (ISCO)

# Bardeen's closed form, x = ±1 and a ≠ 0
function kerr_equatorial_isco(a::Real, x::Real)
    # Use real cube roots to avoid complex branches for negative arguments
    Z1 = 1 + cbrt(1 - a^2) * (cbrt(1 + a) + cbrt(1 - a))
    Z2 = sqrt(3*a^2 + Z1^2)

    # ((3 - Z1)*(3 + Z1 + 2 Z2) / (a x)^2) is non-negative for physical a,x
    denom = (a * x)^2
    inner = ((3 - Z1) * (3 + Z1 + 2*Z2)) / denom

    return 3 + Z2 - (x * a) * sqrt(inner)
end

"""
    kerr_geo_isco(a, x)

Return the ISCO radius: 6 for `a = 0`, and for `a ≠ 0` the equatorial ISCO with `x = ±1`
(prograde when `a x > 0`). For inclined orbits the innermost stable spherical orbit is
`kerr_geo_isso(a, x)`.
"""
function kerr_geo_isco(a::Real, x::Real)
    T = _float_type(a, x)
    if isapprox(a, 0; atol=_tol(T, 1e-12))
        return T(6)
    elseif isapprox(abs(x), 1; atol=_tol(T, 1e-12))
        return kerr_equatorial_isco(a, x)
    else
        throw(DomainError(x, "kerr_geo_isco is defined for a = 0 and for equatorial " *
            "orbits (x = ±1); for inclined orbits use kerr_geo_isso(a, x)."))
    end
end

# Photon Sphere

# x = ±1: prograde for a x > 0
kerr_equatorial_photon_sphere_radius(a::Real, x::Real) =
    2 * (1 + cos(_float_type(a, x)(2) / 3 * acos(-a * sign(x))))

function kerr_polar_photon_sphere_radius(a::Real, x::Real)
    T = _float_type(a, x)
    arg_denom = (1 - (a^2) / 3)
    @assert arg_denom > 0 "Polar formula domain violation (|a| too large for this closed form)."
    inside = (1 - a^2) / (arg_denom^(T(3) / 2))
    inside_clamped = clamp(inside, -1, 1)
    return 1 + 2 * sqrt(arg_denom) * cos(T(1) / 3 * acos(inside_clamped))
end

# |a| = 1 (the caller's test); (a, x) → (−a, −x) is a symmetry
function kerr_extremal_photon_sphere_radius(a::Real, x::Real)
    a < 0 && return kerr_extremal_photon_sphere_radius(-a, -x)
    T = _float_type(a, x)
    return x < sqrt(T(3)) - 1 ? 1 + sqrt(T(2)) * sqrt(1 - x) - x : one(T)
end

function kerr_geo_photon_sphere_radius_numeric(a::Real, x0::Real)
    @assert abs(x0) <= 1.0 "Inclination x must satisfy |x| ≤ 1"
    @assert 0 < abs(a) < 1 "Numeric photon radius requires 0 < |a| < 1."
    # Eliminate the null constants using K = (Phi/x - a*x)^2 and R = R' = 0.
    # With t = r - 1, Delta = (t - delta)(t + delta); no division by a or r - 1
    # is needed. The physical root is bracketed by the outer horizon and r = 4.
    delta2 = (1 - a) * (1 + a)
    delta = sqrt(delta2)
    f(t) = t^3 + (a^2 * (1 + x0^2) - 3) * t - 2 * delta2 +
        2 * a * x0 * (1 + t) * sqrt((t - delta) * (t + delta))
    return 1 + find_zero(f, (delta, 3one(delta)), Bisection())
end

function kerr_geo_photon_sphere_radius(a::Real, x::Real)
    @assert abs(x) <= 1.0 "Inclination parameter x must satisfy |x| ≤ 1"

    T = _float_type(a, x)
    # Schwarzschild case
    if isapprox(a, 0; atol=_tol(T, 1e-12))
        return T(3)
    end

    # Extremal analytic
    if isapprox(abs(a), 1; atol=_tol(T, 1e-12))
        return kerr_extremal_photon_sphere_radius(a, x)
    end

    # Equatorial analytic
    if isapprox(abs(x), 1; atol=_tol(T, 1e-12))
        return kerr_equatorial_photon_sphere_radius(a, x)
    end

    # Polar analytic (x ~ 0)
    if isapprox(abs(x), 0; atol=_tol(T, 1e-14))
        return kerr_polar_photon_sphere_radius(a, x)
    end

    # Generic numeric case
    return kerr_geo_photon_sphere_radius_numeric(a, x)
end

# Separatrix, IBSO and ISSO

"""
    kerr_geo_separatrix(a, e, x)

The separatrix `ps(a, e, x)`: the semi-latus rectum at which the pericentre p/(1 + e) is a
double root of the radial potential. For e < 1 it is the smallest p of a stable orbit with
eccentricity `e` and inclination `x` (below it the orbit plunges; on it the orbit is Critical,
the ISCO or ISSO for e = 0 and a homoclinic orbit for 0 < e < 1); for e ≥ 1 it separates
scattering from capture (2 r_IBSO for e = 1). It is found by bisection on the sign of the
deflated radial polynomial at the pericentre, with the constants of motion of each trial
orbit, between the extremal prograde and retrograde equatorial values 1 + e and
5 + e + 4√(1 + e). When no timelike orbit of this eccentricity and inclination has a double
root at its pericentre (in Schwarzschild spacetime e > 3, where E² < 0 at p = 6 + 2e), a
`DomainError` is raised.
"""
function kerr_geo_separatrix(a::Real, e::Real, x::Real)
    (abs(a) <= 1 && e >= 0 && abs(x) <= 1) || throw(DomainError((a, e, x),
        "kerr_geo_separatrix needs |a| ≤ 1, e ≥ 0 and |x| ≤ 1."))
    lo, hi = float(1 + e), float(5 + e + 4 * sqrt(1 + e))
    # at |a| = 1 the prograde equatorial family ends on the horizon, p = 1 + e, where its constants
    # are degenerate (r = 1 is a double root of Δ and of R) and the residual is rounding noise
    abs(a) == 1 && x == sign(a) && return lo
    # true inside the separatrix, false outside, missing where no timelike orbit has these roots
    inside(p) = try
        _apex_separatrix_residual(a, p, e, x) > 0
    catch err
        err isa DomainError || rethrow()
        missing
    end
    # p = 1 + e (pericentre on the horizon) lies below every separatrix; it is never evaluated,
    # since at |a| = 1 it is the degenerate horizon-root family
    lo_state = missing
    # for e > 1 the family of timelike orbits can begin above the extremal retrograde value
    while inside(hi) !== false
        hi *= 2
        isfinite(hi) || throw(DomainError((a, e, x), "No stable or scattering orbit found."))
    end
    p0 = lo
    while true
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        state = inside(mid)
        if state === false
            hi = mid
        else
            lo, lo_state = mid, state
        end
    end
    lo_state === true && return hi
    # The sign change sits at the lower end of the family of timelike orbits. Either the family ends
    # on a marginally stable orbit at the horizon limit p = 1 + e (|a| = 1, prograde: the constants
    # tend to finite values, and rounding hides the residual within the last ~1e-8 of p) or it
    # ends where E² ∝ 1/(p − p_b) diverges (in Schwarzschild spacetime p_b = 3 + e² for e > 3, the
    # null limit; every timelike orbit with e ≥ 3 scatters). At twice the distance from the lower
    # end the energy is unchanged in the first case and smaller by √2 in the second.
    E_near = _apex_constants(a, hi, e, x).E
    E_far = try
        _apex_constants(a, p0 + 2 * (hi - p0), e, x).E
    catch err
        err isa DomainError || rethrow()
        E_near                  # still within the rounding sliver above the horizon limit
    end
    E_near < 2^(one(E_near) / 4) * E_far && return p0
    throw(DomainError((a, e, x),
        "No timelike orbit with this eccentricity and inclination has a double root at its pericentre."))
end

"""
    kerr_geo_ibso(a, x)

The radius of the innermost bound spherical orbit of inclination `x`: the unstable spherical
orbit with E = 1, also called the marginally bound orbit; half the e = 1 separatrix.
"""
kerr_geo_ibso(a::Real, x::Real) = kerr_geo_separatrix(a, one(_float_type(a, x)), x) / 2

# Innermost stable spherical orbit (ISSO)

"""
    kerr_geo_isso(a, x)

Return the ISSO radius: the innermost stable spherical orbit of inclination `x`, i.e. the
`e = 0` separatrix (the equatorial ISCO for `x = ±1`).
"""
kerr_geo_isso(a::Real, x::Real) = kerr_geo_separatrix(a, zero(_float_type(a, x)), x)

"""
    kerr_geo_orbit_type_metadata(a, p, e, x)

Class of the orbit with APEX-like parameters `(a, p, e, x)`. `family` is the display name of
its class ("Stable", "Critical", "Plunge", "Capture", "Scatter"; "Critical" for the unstable
circular orbits, the ISCO/ISSO and the orbits on the separatrix), `outcome` the class
Symbol, `labels = [family, shape, inclination]`. `stability` ("Stable", "MarginallyStable",
"Unstable", "NotApplicable") and `energy_regime` ("Elliptic", "Parabolic", "Hyperbolic") are
separate fields. Parameters within `separatrix_tolerance` of the separatrix are put on it.
For `p ≤ 0` the orbit is not classified: `family = "NotClassified"`, `outcome = :not_classified`,
`stability = "Unknown"` and `labels = [family]`.
"""
function kerr_geo_orbit_type_metadata(a::Real, p::Real, e::Real, x::Real)
    T = _float_type(a, p, e, x)
    iszero_tol(v) = isapprox(v, 0; atol=_tol(T, 1e-12))
    circular = iszero_tol(e)
    inclination = iszero_tol(abs(x) - 1) ? "Equatorial" : "Inclined"
    shape = circular ? "Circular" : e < 1 ? "Eccentric" : iszero_tol(e - 1) ? "Parabolic" : "Hyperbolic"
    # e ≥ 1 without a separatrix (no timelike orbit of this eccentricity has a double root at
    # its pericentre): every orbit scatters
    separatrix_p = circular ? kerr_geo_isso(a, x) :
        try kerr_geo_separatrix(a, e, x) catch err
            (err isa DomainError && e >= 1) || rethrow()
            T(NaN)
        end
    photon_p = circular ? kerr_geo_photon_sphere_radius(a, x) : T(NaN)
    ibso_p = circular ? kerr_geo_ibso(a, x) : T(NaN)
    isso_p = circular ? kerr_geo_isso(a, x) : T(NaN)
    tolerance = _tol(T, 1e-12)
    separatrix_tolerance = _tol(T, 1e-15)
    at_separatrix = abs(p - separatrix_p) <= separatrix_tolerance
    p_effective = at_separatrix ? separatrix_p : p
    on_separatrix = isapprox(p_effective, separatrix_p; atol=tolerance)
    # E < 1, = 1, > 1: from e, and for circular orbits from p against the IBSO
    energy_regime = circular ?
        (isapprox(p_effective, ibso_p; atol=tolerance) ? "Parabolic" :
         p_effective < ibso_p ? "Hyperbolic" : "Elliptic") :
        e < 1 ? "Elliptic" : iszero_tol(e - 1) ? "Parabolic" : "Hyperbolic"
    family, stability = if p_effective <= 0
        "NotClassified", "Unknown"
    elseif circular
        p_effective <= photon_p ? ("Plunge", "Unstable") :
        on_separatrix ? ("Critical", "MarginallyStable") :
        p_effective < isso_p ? ("Critical", "Unstable") : ("Stable", "Stable")
    elseif e < 1
        # the eccentric separatrix orbit is homoclinic to an unstable spherical orbit
        on_separatrix ? ("Critical", "Unstable") :
        p_effective > separatrix_p ? ("Stable", "Stable") : ("Plunge", "Unstable")
    else
        on_separatrix ? ("Critical", "Unstable") :
        p_effective > separatrix_p || isnan(separatrix_p) ? ("Scatter", "NotApplicable") :
        ("Capture", "Unstable")
    end
    known = family != "NotClassified"
    return (
        labels=known ? [family, shape, inclination] : [family],
        family=family,
        outcome=known ? Symbol(lowercase(family)) : :not_classified,
        shape=shape,
        inclination=inclination,
        energy_regime=energy_regime,
        stability=stability,
        input_p=p,
        effective_p=p_effective,
        separatrix_p=separatrix_p,
        at_separatrix=at_separatrix,
        photon_p=photon_p,
        ibso_p=ibso_p,
        isso_p=isso_p,
        tolerance=tolerance,
        separatrix_tolerance=separatrix_tolerance,
    )
end

"""
    kerr_geo_orbit_type(a, p, e, x)

Return the orbit-type labels `[family, shape, inclination]` of
`kerr_geo_orbit_type_metadata(a, p, e, x)`.
"""
function kerr_geo_orbit_type(a::Real, p::Real, e::Real, x::Real)
    return kerr_geo_orbit_type_metadata(a, p, e, x).labels
end
