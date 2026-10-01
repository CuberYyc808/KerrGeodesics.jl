# APEX reference API: constants of motion (E, Lz, Q) from (a, p, e, x) through the solver of
# models/ApexConstants.jl; the Schwarzschild scatter quantities (v∞, b, ψ) are the closed forms of
# the Mathematica KerrGeodesics package.

"""
    kerr_geo_energy(a, p, e, x)

Energy E of the orbit with APEX parameters `(a, p, e, x)`: r = p/(1 ± e) are roots of the
radial potential and x = cos of the inclination, Q = (1 − x²)(a²(1 − E²) + Lz²/x²).
"""
kerr_geo_energy(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).E

"""
    kerr_geo_angular_momentum(a, p, e, x)

Axial angular momentum Lz of the orbit with APEX parameters `(a, p, e, x)`; see
[`kerr_geo_energy`](@ref).
"""
kerr_geo_angular_momentum(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).Lz

"""
    kerr_geo_carter_constant(a, p, e, x)

Carter constant Q = (1 − x²)(a²(1 − E²) + Lz²/x²) of the orbit with APEX parameters
`(a, p, e, x)`; see [`kerr_geo_energy`](@ref).
"""
kerr_geo_carter_constant(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).Q

# Scatter quantities of a Schwarzschild (a = 0) orbit with e > 1.
"""
    kerr_geo_velocity_at_infinity(a, p, e, x)

Speed at infinity v∞ = √(E² − 1)/E of a Schwarzschild (a = 0) orbit with e > 1.
"""
function kerr_geo_velocity_at_infinity(a::Real, p::Real, e::Real, x::Real)
    iszero(a) || throw(ArgumentError("kerr_geo_velocity_at_infinity: the formula holds for a = 0 only"))
    e > 1 || throw(ArgumentError("kerr_geo_velocity_at_infinity: e = $e; an orbit reaching infinity with v∞ > 0 has e > 1"))
    En = kerr_geo_energy(a, p, e, x)
    return sqrt(En^2 - 1) / En
end

"""
    kerr_geo_impact_parameter(a, p, e, x)

Impact parameter b = Lz/√(E² − 1) of a Schwarzschild (a = 0) orbit with e > 1.
"""
function kerr_geo_impact_parameter(a::Real, p::Real, e::Real, x::Real)
    iszero(a) || throw(ArgumentError("kerr_geo_impact_parameter: the formula holds for a = 0 only"))
    e > 1 || throw(ArgumentError("kerr_geo_impact_parameter: e = $e; an orbit reaching infinity with v∞ > 0 has e > 1"))
    Lz = kerr_geo_angular_momentum(a, p, e, x)
    En = kerr_geo_energy(a, p, e, x)
    return Lz / sqrt(En^2 - 1)
end

"""
    kerr_geo_deflection_angle(a, p, e, x=1.0)

Deflection angle ψ = Δφ − π of a Schwarzschild (a = 0) orbit from infinity through the
pericentre r = p/(1 + e) back to infinity (e ≥ 1; e = 1 is the parabolic orbit). With
u = 1/r the orbit equation is (du/dφ)² = 2(u₁ − u)(u − u₂)(u₃ − u), u₁ = (1 + e)/p,
u₂ = (1 − e)/p, u₃ = 1/2 − 2/p; u = u₂ + (u₁ − u₂) sin²χ gives

    ψ = 4√(p/Δ) [K(k) − F(χ₀, k)] − π,   Δ = p − 6 + 2e,   k = 4e/Δ,   cos 2χ₀ = 1/e,

with K, F in parameter convention. The pericentre is a turning point reached from infinity
only for u₁ < u₃, i.e. p > 6 + 2e (k < 1); for p ≤ 6 + 2e the orbit is captured and there is
no deflection angle. A timelike orbit also needs p > 3 + e² (L² = p²/(p − 3 − e²)), which is
the stronger condition for e > 3. For p → ∞, ψ = 2 asin(1/e) + [6 acos(−1/e) + 2√(e² − 1)]/p + O(1/p²).
"""
function kerr_geo_deflection_angle(a::Real, p::Real, e::Real, x::Real=1.0)
    iszero(a) || throw(ArgumentError("kerr_geo_deflection_angle: the closed form holds for a = 0 only"))
    e >= 1 || throw(ArgumentError("kerr_geo_deflection_angle: e = $e < 1 is a bound orbit"))
    p > 6 + 2e || throw(ArgumentError("kerr_geo_deflection_angle: p = $p ≤ 6 + 2e for e = $e; " *
        "the orbit is captured and has no deflection angle"))
    p > 3 + e^2 || throw(ArgumentError("kerr_geo_deflection_angle: p = $p ≤ 3 + e² for e = $e; " *
        "no timelike orbit has these (p, e), since L² = p²/(p − 3 − e²)"))
    Δ = p - 6 + 2e
    k = 4e / Δ
    χ0 = acos(1 / e) / 2
    return 4 * sqrt(p / Δ) * (Elliptic.K(k) - Elliptic.F(χ0, k)) - π
end

"""
    kerr_geo_constants_of_motion(a, p, e, x)

The constants of motion of the orbit with APEX parameters `(a, p, e, x)`, as
`Dict("E" => E, "Lz" => Lz, "Q" => Q)`; r = p/(1 ± e) are roots of the radial potential and
x = cos of the inclination, Q = (1 − x²)(a²(1 − E²) + Lz²/x²). A Schwarzschild (a = 0) scatter
orbit (e > 1 and p > 6 + 2e, so it returns to infinity) also carries the speed at infinity
"v∞", the impact parameter "b" and the deflection angle "ψ".

The parameters must admit real timelike constants, otherwise a `DomainError` is raised. A root
of the separatrix polynomial need not satisfy this at large eccentricity: in Schwarzschild
spacetime `p = 6 + 2e` with `e = 5` gives E² = −1/2. Scattering with `e = 5` exists at larger
`p`; in Schwarzschild spacetime it requires `p > 3 + e²` as well as `p > 6 + 2e`.
"""
function kerr_geo_constants_of_motion(a::Real, p::Real, e::Real, x::Real)
    c = _apex_constants(a, p, e, x)
    constants = Dict("E" => c.E, "Lz" => c.Lz, "Q" => c.Q)
    if iszero(a) && e > 1 && p > 6 + 2e
        constants["v∞"] = kerr_geo_velocity_at_infinity(a, p, e, x)
        constants["b"] = kerr_geo_impact_parameter(a, p, e, x)
        constants["ψ"] = kerr_geo_deflection_angle(a, p, e, x)
    end
    return constants
end
