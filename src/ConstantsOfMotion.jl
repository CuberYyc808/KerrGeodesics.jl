module ConstantsOfMotion

# This module translates the Mathematica ConstantsOfMotion subpackage
# into Julia. It implements Schwarzschild (a=0) special cases and the
# Kerr cases (equatorial, polar, spherical, marginally bound, generic).

using Elliptic
export kerr_geo_energy,
        kerr_geo_angular_momentum,
        kerr_geo_carter_constant,
        kerr_geo_constants_of_motion,
        kerr_geo_velocity_at_infinity,
        kerr_geo_impact_parameter,
        kerr_geo_scattering_angle

const ATOL = 1e-12

#########################
# Helper utilities
#########################

_is_one_sq(x::Real) = isapprox(x^2, 1.0; atol=ATOL)
_is_zero(x::Real) = isapprox(x, 0.0; atol=ATOL)

# Schwarzschild special-case wrappers (a == 0)
# We'll implement these as helper functions named schwarz_* for clarity,
# then the main kerr_* functions will dispatch to them when a == 0.

#########################
# Schwarzschild formulas
#########################

"""
    schwarz_geo_energy(p, e, x)

Schwarzschild orbital energy:
- Circular (e=0): E = (p-2)/sqrt(p*(p-3))
- Eccentric: E = sqrt(((p-2)^2 - 4 e^2) / (p*(p-3-e^2)))
"""
function schwarz_geo_energy(p::Real, e::Real, x::Real)
    if _is_one_sq(x) && _is_zero(e)  # circular equatorial
        return (p - 2) / sqrt(p * (p - 3))
    else
        # generic eccentric Schwarzschild bound formula (from your MMA)
        return sqrt(((-4*e^2 + ( -2 + p )^2) ) / (p * ( -3 - e^2 + p )))
    end
end

"""
    schwarz_geo_angular_momentum(p, e, x)

Schwarzschild L_z:
- Circular (e=0): L = p * x / sqrt(p-3)
- Eccentric: L = p * x / sqrt(p - 3 - e^2)
"""
function schwarz_geo_angular_momentum(p::Real, e::Real, x::Real)
    if _is_zero(e)
        return p * x / sqrt(p - 3)
    else
        return p * x / sqrt(p - 3 - e^2)
    end
end

"""
    schwarz_geo_carter_constant(p, e, x)

Schwarzschild Carter constant:
- circular: Q = - p^2 (1 - x^2)/( -3 + p )
- eccentric: Q = p^2 (1 - x^2)/( 3 + e^2 - p )
Note: the sign/representation matches the Mathematica form.
"""
function schwarz_geo_carter_constant(p::Real, e::Real, x::Real)
    if _is_zero(e)
        return - (p^2 * (1 - x^2)) / ( -3 + p )
    else
        return (p^2 * (1 - x^2)) / ( 3 + e^2 - p )
    end
end

# Scatter helpers (Schwarzschild, e > 1)
"""
    kerr_geo_velocity_at_infinity(a, p, e, x)

v_infty = sqrt(E^2 - 1) / E  (defined only for scatter e > 1).
If called with a != 0, it will dispatch to Kerr formulas; here we just
implement Schwarzschild mapping (a==0).
"""

# Energy at infinity for hyperbolic (e > 1)
function kerr_geo_velocity_at_infinity(a::Real, p::Real, e::Real, x::Real)
    @assert isapprox(a, 0.0; atol=1e-12) "This function is for Schwarzschild (a=0)"
    @assert e > 1 "Hyperbolic motion requires e > 1"
    En = kerr_geo_energy(a, p, e, x)   # You need to define kerr_geo_energy for Schwarzschild
    return sqrt(En^2 - 1) / En
end

# Impact parameter for hyperbolic motion
function kerr_geo_impact_parameter(a::Real, p::Real, e::Real, x::Real)
    @assert isapprox(a, 0.0; atol=1e-12) "This function is for Schwarzschild (a=0)"
    @assert e > 1 "Hyperbolic motion requires e > 1"
    Lz = kerr_geo_angular_momentum(a, p, e, x)  # You need to define kerr_geo_angular_momentum
    En = kerr_geo_energy(a, p, e, x)
    return Lz / sqrt(En^2 - 1)
end

# Schwarzschild hyperbolic scattering angle
function kerr_geo_scattering_angle(a::Real, p::Real, e::Real, x::Real=1.0)
    @assert isapprox(a, 0.0; atol=1e-12) "This formula is for Schwarzschild (a=0)"
    @assert e >= 1 "Scattering angle defined for e >= 1"

    # Check that orbit is physically allowed: denominator must be positive
    Δ = p - 6 - 2*e
    if Δ <= 0
        return error("Orbit parameters invalid: p=$p is too small for e=$e, particle falls into the BH.")
    end

    # Modulus for elliptic integral
    m = (4*e) / Δ
    
    # Amplitude
    phi = acos(-1 / e)
    
    # Scattering angle
    θ = -π + (4 * sqrt(p) * Elliptic.F(phi, m)) / sqrt(Δ)
    
    # Check if result is real
    if !isreal(θ)
        return error("p=$p is too small for e=$e as a scattering orbit")
    end

    return real(θ)
end


"""
    schwarz_geo_constants(p, e, x)

Return NamedTuple of Schwarzschild constants (E, Lz, Q) or with scattering extra fields if e>=1.
"""
function schwarz_geo_constants(p::Real, e::Real, x::Real)
    E = schwarz_geo_energy(p, e, x)
    Lz = schwarz_geo_angular_momentum(p, e, x)
    Q = schwarz_geo_carter_constant(p, e, x)
    if e >= 1
        vinf = kerr_geo_velocity_at_infinity(0.0, p, e, x)
        b = kerr_geo_impact_parameter(0.0, p, e, x)
        ψ = kerr_geo_scattering_angle(0.0, p, e, x)
        return Dict("E"=>E, "Lz"=>Lz, "Q"=>Q, "v∞"=>vinf, "b"=>b, "ψ"=>ψ)
    else
        return Dict("E"=>E, "Lz"=>Lz, "Q"=>Q)
    end
end

#########################
# Kerr formulas
#########################

# Carter constant is zero for equatorial orbits (x^2 == 1)
# We'll use isapprox for testing x^2==1.

# Equatorial circular (e=0, x^2==1)
function _kerr_energy_equatorial_circular(a::Real, p::Real, x::Real)
    return ((-2 + p) * sqrt(p) + a / x) / sqrt( (2*a/x) * p^(3/2) + ( -3 + p ) * p^2 )
end

function _kerr_L_equatorial_circular(a::Real, p::Real, x::Real)
    num = ((a^2 + p^2) * x - 2*a*sqrt(p))
    den = p^(3/4) * sqrt( x^2 * (-3 + p) * sqrt(p) + 2*a*x )
    return num / den
end

# Equatorial eccentric (from Glampedakis & Kennefick simplification)
function _kerr_energy_equatorial_eccentric(a::Real, p::Real, e::Real, x::Real)
    A = a^2 * (1 + 3*e^2 + p)
    B = p * (-3 - e^2 + p)
    inner_sqrt = sqrt( ( a^6*( -1 + e^2 )^2 + a^2*( -4*e^2 + ( -2 + p )^2 ) * p^2 + 2*a^4*p*( -2 + p + e^2*(2 + p) ) ) / (p^3 * x^2) )
    numer = ( -1 + e^2 ) * ( A + B - 2*x*inner_sqrt )
    denom = -4*a^2*( -1 + e^2 )^2 + (3 + e^2 - p)^2 * p
    # outer expression as in MMA:
    return sqrt( 1 - ( (1 - e^2) * ( 1 + ( (-1 + e^2) * ( A + p*( -3 - e^2 + p - 2*x*inner_sqrt ) ) ) / denom ) ) / p )
end

function _kerr_L_equatorial_eccentric(a::Real, p::Real, e::Real, x::Real)
    # This matches the MMA expression:
    inner_sqrt = sqrt( ( a^6*( -1 + e^2 )^2 + a^2*( -4*e^2 + ( -2 + p )^2 ) * p^2 + 2*a^4*p*( -2 + p + e^2*(2 + p) ) ) / ( p^3 * x^2 ) )
    numerator1 = a^2*(1 + 3*e^2 + p) + p*( -3 - e^2 + p - 2*x*inner_sqrt )
    denom1 = ( -4*a^2*( -1 + e^2 )^2 + (3 + e^2 - p)^2 * p ) * x^2
    part1 = p * x * sqrt( numerator1 / denom1 )
    # the + a * sqrt[...] term (same radial energy-like sqrt)
    part2 = a * sqrt( 1 - ( (1 - e^2) * ( 1 + ( (-1 + e^2) * ( a^2*(1 + 3*e^2 + p) + p*( -3 - e^2 + p - 2*x*inner_sqrt ) ) ) / ( -4*a^2*( -1 + e^2 )^2 + (3 + e^2 - p)^2 * p ) ) ) / p )
    return part1 + part2
end

# Polar orbits (x == 0): Lz = 0
function _is_polar(x::Real)
    return isapprox(x, 0.0; atol=ATOL)
end

#########################
# Spherical orbits (x arbitrary, e = 0, "spherical" means r=const different from equatorial)
# There are two variants: a special-case (both e=0 and x=0) with simple forms and
# a general spherical (e=0, arbitrary x) with long expressions.
#########################

# Spherical, special case x=0, e=0 (from Stoghianidis & Tsoubelis)
function _kerr_energy_spherical_polar(a::Real, p::Real)
    num = p * (a^2 - 2*p + p^2)^2
    den = (a^2 + p^2) * (a^2 + a^2*p - 3*p^2 + p^3)
    return sqrt( num / den )
end

function _kerr_Q_spherical_polar(a::Real, p::Real)
    return ( p^2 * ( a^4 + 2*a^2*( -2 + p ) * p + p^4 ) ) / ( (a^2 + p^2) * ( (-3 + p)*p^2 + a^2*(1 + p) ) )
end

# General spherical orbits (e=0, arbitrary x)
function _kerr_energy_spherical_general(a::Real, p::Real, x::Real)
    # To avoid an unreadably long single expression, compute subterms:
    T0 = (-3 + p) * (-2 + p)^2 * p^5
    T1 = -2 * a^5 * x * (x^2 - 1) * sqrt(p^3 + a^2 * p * ( -1 + x^2 ))
    T2 = a^4 * p^2 * ( -1 + x^2 ) * ( 4 - 5*p*( -1 + x^2 ) + 3*p^2*( -1 + x^2 ) )
    T3 = - a^6 * ( -1 + x^2 )^2 * ( x^2 + p^2*( -1 + x^2 ) - p * (1 + 2*x^2) )
    T4 = a^2 * p^3 * ( 4 - 4*x^2 + p*(12 - 7*x^2) - 3*p^3*( -1 + x^2 ) + p^2*( -13 + 10*x^2 ) )
    T5 = a * ( -2 * p^(9/2) * x * sqrt( p^2 + a^2 * ( -1 + x^2 ) ) + 4 * p^3 * x * sqrt( p^3 + a^2 * p * ( -1 + x^2 ) ) )
    T6 = 2 * a^3 * ( 2*p*x*( -1 + x^2 ) * sqrt( p^3 + a^2 * p * ( -1 + x^2 ) ) - x^3 * sqrt( p^7 + a^2 * p^5 * ( -1 + x^2 ) ) )

    numerator = T0 + T1 + T2 + T3 + T4 + T5 + T6

    denom1 = ( p^2 - a^2 * ( -1 + x^2 ) )
    denom2 = ( (-3 + p)^2 * p^4 - 2*a^2 * p^2 * (3 + 2*p - 3*x^2 + p^2*( -1 + x^2 )) + a^4 * ( -1 + x^2 ) * ( -1 + x^2 + p^2*( -1 + x^2 ) - 2*p*(1 + x^2) ) )

    return sqrt( numerator / ( denom1 * denom2 ) )
end

# Angular momentum from spherical solution (uses En optionally)
function _kerr_L_from_spherical(a::Real, p::Real, x::Real; En1=nothing)
    En = En1 === nothing ? kerr_geo_energy(a, p, 0.0, x) : En1
    g = 2 * a * p
    d = ( a^2 + ( -2 + p ) * p ) * ( p^2 - a^2 * ( -1 + x^2 ) )
    h = ( (-2 + p) * p - a^2 * ( -1 + x^2 ) ) / ( x^2 )
    f = p^4 + a^2 * ( p*(2 + p) - ( a^2 + ( -2 + p ) * p ) * ( -1 + x^2 ) )

    val = (-En * g + x * sqrt( ( -d*h + En^2*( g^2 + f*h ) ) / x^2 )) / h
    return val
end

#########################
# Marginally bound (e == 1)
#########################
function _kerr_energy_marginally_bound(a::Real, p::Real, x::Real)
    return 1.0
end

function _kerr_L_marginally_bound(a::Real, p::Real, x::Real; En1=nothing)
    En = En1 === nothing ? kerr_geo_energy(a, p, 1.0, x) : En1
    ρ2 = p / (1 + 1.0)  # since e == 1
    numerator = x^2 * ( -2*a*ρ2 + sqrt(2)/x * sqrt( ρ2 * ( a^2 + ( -2 + ρ2 ) * ρ2 ) * ( a^2*(1 - x^2) + ρ2^2 ) ) )
    denom = a^2 * (1 - x^2) + ( -2 + ρ2 ) * ρ2
    return numerator / denom
end

#########################
# Generic case (most general formulas for Kerr)
#########################

# Generic Kerr energy:
function kerr_geo_energy(a::Real, p::Real, e::Real, x::Real)
    # If a == 0, dispatch to Schwarzschild versions
    if isapprox(a, 0.0; atol=ATOL)
        return schwarz_geo_energy(p, e, x)
    end

    # Negative-a mapping: use symmetry a<0 -> transform to (-a, -x)
    if a < 0
        return kerr_geo_energy(-a, p, e, -x)
    end

    # Equatorial special cases
    if _is_one_sq(x)
        if _is_zero(e)
            return _kerr_energy_equatorial_circular(a, p, x)
        else
            return _kerr_energy_equatorial_eccentric(a, p, e, x)
        end
    end

    # Polar special-case: x == 0
    if _is_polar(x)
        if _is_zero(e)
            # spherical polar special case simplified
            return _kerr_energy_spherical_polar(a, p)
        else
            # eccentric polar case (Warburton / Schmidt derived)
            numer = - ( p * ( a^4*( -1 + e^2 )^2 + ( -4*e^2 + ( -2 + p )^2 ) * p^2 + 2*a^2*p*( -2 + p + e^2*(2 + p) ) ) )
            denom = a^4*( -1 + e^2 )^2 * ( -1 + e^2 - p ) + (3 + e^2 - p) * p^4 - 2*a^2 * p^2 * ( -1 - e^4 + p + e^2*(2 + p) )
            return sqrt( numer / denom )
        end
    end

    # Spherical orbits with e == 0 and arbitrary x
    if _is_zero(e)
        return _kerr_energy_spherical_general(a, p, x)
    end

    # Marginally bound e == 1
    if isapprox(e, 1.0; atol=ATOL)
        return _kerr_energy_marginally_bound(a, p, x)
    end

    # Generic bound/eccentric case (most general formula from MMA)
    # Implementing the module body that computes r1,r2,zm and the f,g,h,d functions,
    # then constructs the coefficients and returns the final sqrt expression.
    r1 = p / (1 - e)
    r2 = p / (1 + e)
    zm = sqrt(1 - x^2)

    Delta(r) = r^2 - 2*r + a^2

    f(r) = r^4 + a^2 * ( r*(r + 2) + zm^2 * Delta(r) )
    g(r) = 2 * a * r
    h(r) = r*(r - 2) + (zm^2 / (1 - zm^2)) * Delta(r)
    d(r) = ( r^2 + a^2 * zm^2 ) * Delta(r)

    κ = d(r1) * h(r2) - h(r1) * d(r2)
    ε = d(r1) * g(r2) - g(r1) * d(r2)
    ρ = f(r1) * h(r2) - h(r1) * f(r2)
    η = f(r1) * g(r2) - g(r1) * f(r2)
    σ = g(r1) * h(r2) - h(r1) * g(r2)

    # inner sqrt term (note x could be small; caller should avoid x==0 here)
    inner = σ * ( σ * ε^2 + ρ * ε * κ - η * κ^2 )
    # main numerator & denominator
    numer = κ * ρ + 2 * ε * σ - 2 * x * sqrt( inner / x^2 )
    denom = ρ^2 + 4 * η * σ

    return sqrt( numer / denom )
end

# Generic angular momentum L_z
function kerr_geo_angular_momentum(a::Real, p::Real, e::Real, x::Real; En1::Union{Nothing,Real}=nothing)
    # Schwarzschild dispatch
    if isapprox(a, 0.0; atol=ATOL)
        return schwarz_geo_angular_momentum(p, e, x)
    end

    # negative-a
    if a < 0
        return -kerr_geo_angular_momentum(-a, p, e, -x; En1=En1)
    end

    # equatorial
    if _is_one_sq(x)
        if _is_zero(e)
            return _kerr_L_equatorial_circular(a, p, x)
        else
            return _kerr_L_equatorial_eccentric(a, p, e, x)
        end
    end

    # polar
    if _is_polar(x)
        return 0.0
    end

    # spherical e==0 but arbitrary x: use spherical L compute (requires En)
    if _is_zero(e)
        En = En1 === nothing ? kerr_geo_energy(a, p, e, x) : En1
        return _kerr_L_from_spherical(a, p, x; En1=En)
    end

    # marginally bound e==1
    if isapprox(e, 1.0; atol=ATOL)
        En = En1 === nothing ? kerr_geo_energy(a, p, e, x) : En1
        return _kerr_L_marginally_bound(a, p, x; En1=En)
    end

    # generic case
    En = En1 === nothing ? kerr_geo_energy(a, p, e, x) : En1
    r1 = p / (1 - e)
    zm = sqrt(1 - x^2)

    Delta(r) = r^2 - 2*r + a^2
    f(r) = r^4 + a^2 * ( r*(r + 2) + zm^2 * Delta(r) )
    g(r) = 2 * a * r
    h(r) = r*(r - 2) + (zm^2 / (1 - zm^2)) * Delta(r)
    d(r) = ( r^2 + a^2 * zm^2 ) * Delta(r)

    return (-En * g(r1) + x * sqrt( ( -d(r1) * h(r1) + En^2 * ( g(r1)^2 + f(r1) * h(r1) ) ) / x^2 )) / h(r1)
end

# Carter constant
function kerr_geo_carter_constant(a::Real, p::Real, e::Real, x::Real; En1::Union{Nothing,Real}=nothing, L1::Union{Nothing,Real}=nothing)
    # Schwarzschild dispatch
    if isapprox(a, 0.0; atol=ATOL)
        return schwarz_geo_carter_constant(p, e, x)
    end

    if isapprox(x, 0.0; atol=1e-12)
        numerator = -p^2 * (a^4*(1 - e^2)^2 + p^4 + 2*a^2*p*(-2 + p + e^2*(2+p)))
        denominator = a^4*(1 - e^2)^2*( -1 + e^2 - p ) + (3 + e^2 - p)*p^4 - 2*a^2*p^2*( -1 - e^4 + p + e^2*(2+p) )
        return numerator / denominator
    end

    En = En1 === nothing ? kerr_geo_energy(a, p, e, x) : En1
    Lz = L1 === nothing ? kerr_geo_angular_momentum(a, p, e, x; En1=En) : L1

    zm = sqrt(1 - x^2)
    return zm^2 * ( a^2 * (1 - En^2) + Lz^2 / (1 - zm^2) )
end

function kerr_geo_constants_of_motion(a::Real, p::Real, e::Real, x::Real)
    if isapprox(a, 0.0; atol=ATOL)
        return schwarz_geo_constants(p, e, x)
    end
    En = kerr_geo_energy(a, p, e, x)
    Lz = kerr_geo_angular_momentum(a, p, e, x; En1=En)
    Q = kerr_geo_carter_constant(a, p, e, x; En1=En, L1=Lz)
    return Dict("E"=>En, "Lz"=>Lz, "Q"=>Q)
end

end
