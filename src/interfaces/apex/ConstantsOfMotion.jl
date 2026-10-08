# APEX reference API: constants of motion (E, Lz, Q) from (a, p, e, x) through the solver of
# models/ApexConstants.jl; the Schwarzschild scatter quantities (v∞, b, ψ) are the closed forms of
# the Mathematica KerrGeodesics package.

"""
    kerr_geo_energy(a, p, e, x)

Energy ``E`` of the orbit with APEX parameters ``(a, p, e, x)``: ``r = p/(1 \\pm e)`` are roots of the
radial potential and ``x`` is the cosine of the inclination,
``Q = (1 - x^2)(a^2(1 - E^2) + L_z^2/x^2)``.
"""
kerr_geo_energy(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).E

"""
    kerr_geo_angular_momentum(a, p, e, x)

Axial angular momentum ``L_z`` of the orbit with APEX parameters ``(a, p, e, x)``; see
[`kerr_geo_energy`](@ref).
"""
kerr_geo_angular_momentum(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).Lz

"""
    kerr_geo_carter_constant(a, p, e, x)

Carter constant ``Q = (1 - x^2)(a^2(1 - E^2) + L_z^2/x^2)`` of the orbit with APEX parameters
``(a, p, e, x)``; see [`kerr_geo_energy`](@ref).
"""
kerr_geo_carter_constant(a::Real, p::Real, e::Real, x::Real) = _apex_constants(a, p, e, x).Q

# Scatter quantities of a Schwarzschild (a = 0) orbit with e > 1.
"""
    kerr_geo_velocity_at_infinity(a, p, e, x)

Speed at infinity ``v_\\infty = \\sqrt{E^2 - 1}/E`` of a Schwarzschild (``a = 0``) orbit with ``e > 1``.
"""
function kerr_geo_velocity_at_infinity(a::Real, p::Real, e::Real, x::Real)
    iszero(a) || throw(ArgumentError("kerr_geo_velocity_at_infinity: the formula holds for a = 0 only"))
    e > 1 || throw(ArgumentError("kerr_geo_velocity_at_infinity: e = $e; an orbit reaching infinity with v∞ > 0 has e > 1"))
    En = kerr_geo_energy(a, p, e, x)
    return sqrt(En^2 - 1) / En
end

"""
    kerr_geo_impact_parameter(a, p, e, x)

Impact parameter ``b = L_z/\\sqrt{E^2 - 1}`` of a Schwarzschild (``a = 0``) orbit with ``e > 1``.
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

Deflection angle ``\\psi = \\Delta\\phi - \\pi`` of a Schwarzschild (``a = 0``) orbit from infinity through the
pericentre ``r = p/(1 + e)`` back to infinity (``e \\geq 1``; ``e = 1`` is the parabolic orbit). With
``u = 1/r`` the orbit equation is ``(du/d\\phi)^2 = 2(u_1 - u)(u - u_2)(u_3 - u)``,
``u_1 = (1 + e)/p``, ``u_2 = (1 - e)/p``, ``u_3 = 1/2 - 2/p``;
``u = u_2 + (u_1 - u_2)\\sin^2\\chi`` gives

```math
\\psi = 4\\sqrt{p/\\Delta} [K(k) - F(\\chi_0, k)] - \\pi, \\qquad
\\Delta = p - 6 + 2e, \\qquad k = 4e/\\Delta, \\qquad \\cos(2\\chi_0) = 1/e.
```

with ``K``, ``F`` in parameter convention. The pericentre is a turning point reached from infinity
only for ``u_1 < u_3``, i.e. ``p > 6 + 2e`` (``k < 1``); for ``p \\leq 6 + 2e`` the orbit is captured and there is
no deflection angle. A timelike orbit also needs ``p > 3 + e^2`` (``L^2 = p^2/(p - 3 - e^2)``), which is
the stronger condition for ``e > 3``. For ``p \\to \\infty``,
``\\psi = 2\\arcsin(1/e) + [6\\arccos(-1/e) + 2\\sqrt{e^2 - 1}]/p + O(1/p^2)``.
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
    return 4 * sqrt(p / Δ) * (_K(k) - _F(χ0, k)) - π
end

"""
    kerr_geo_constants_of_motion(a, p, e, x; precision=nothing)

The constants of motion of the orbit with APEX parameters `(a, p, e, x)`, as
`Dict("E" => E, "Lz" => Lz, "Q" => Q)`; ``r = p/(1 \\pm e)`` are roots of the radial potential and
``x`` is the cosine of the inclination, ``Q = (1 - x^2)(a^2(1 - E^2) + L_z^2/x^2)``.
A Schwarzschild (``a = 0``) scatter orbit (``e > 1`` and ``p > 6 + 2e``, so it returns to infinity)
also carries the speed at infinity `"v∞"`, the impact parameter `"b"` and the deflection angle `"ψ"`.

The parameters must admit real timelike constants, otherwise a `DomainError` is raised. A root
of the separatrix polynomial need not satisfy this at large eccentricity: in Schwarzschild
spacetime ``p = 6 + 2e`` with ``e = 5`` gives ``E^2 = -1/2``. Scattering with ``e = 5`` exists at larger
``p``; in Schwarzschild spacetime it requires ``p > 3 + e^2`` as well as ``p > 6 + 2e``. The constants
are computed in the floating-point type of `(a, p, e, x)`; `precision = p` converts them to
`BigFloat` of `p` bits.
"""
function kerr_geo_constants_of_motion(a::Real, p::Real, e::Real, x::Real; precision=nothing)
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geo_constants_of_motion(BigFloat(a), BigFloat(p), BigFloat(e), BigFloat(x))
    end
    T = _float_type(a, p, e, x)
    return _with_precision(T, _input_precision(a, p, e, x)) do
        _kerr_geo_constants_of_motion(T(a), T(p), T(e), T(x))
    end
end

function _kerr_geo_constants_of_motion(a, p, e, x)
    c = _apex_constants(a, p, e, x)
    constants = Dict("E" => c.E, "Lz" => c.Lz, "Q" => c.Q)
    if iszero(a) && e > 1 && p > 6 + 2e
        constants["v∞"] = kerr_geo_velocity_at_infinity(a, p, e, x)
        constants["b"] = kerr_geo_impact_parameter(a, p, e, x)
        constants["ψ"] = kerr_geo_deflection_angle(a, p, e, x)
    end
    return constants
end
