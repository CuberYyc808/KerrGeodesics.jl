module SpecialOrbits

include("ConstantsOfMotion.jl")
using .ConstantsOfMotion  
using Roots              

export kerr_geo_photon_sphere_radius,
        kerr_geo_isco,
        kerr_geo_ibso,
        kerr_geo_isso,
        kerr_geo_separatrix,
        kerr_geo_orbit_type


# Innermost stable circular orbit (ISCO)

schwarzschild_geo_isco(a::Real, x::Real) = 6.0

function kerr_equatorial_isco(a::Real, x::Real)
    @assert isapprox(abs(x), 1.0; atol=1e-12) "Equatorial ISCO requires x = ±1"
    @assert !isapprox(a, 0.0; atol=1e-12) "For a ≈ 0 use schwarzschild_geo_isco"

    # Use real cube roots to avoid complex branches for negative arguments
    Z1 = 1 + cbrt(1 - a^2) * (cbrt(1 + a) + cbrt(1 - a))
    Z2 = sqrt(3*a^2 + Z1^2)

    # ((3 - Z1)*(3 + Z1 + 2 Z2) / (a x)^2) is non-negative for physical a,x
    denom = (a * x)^2
    inner = ((3 - Z1) * (3 + Z1 + 2*Z2)) / denom

    return 3 + Z2 - (x * a) * sqrt(inner)
end

function kerr_geo_isco(a::Real, x::Real)
    if isapprox(a, 0.0; atol=1e-12)
        return schwarzschild_geo_isco(a, x)
    elseif isapprox(abs(x), 1.0; atol=1e-12)
        return kerr_equatorial_isco(a, x)
    else
        error("kerr_geo_isco: general inclined Kerr ISCO not implemented yet. Implement numeric solver for arbitrary x if required.")
    end
end

# Photon Sphere

schwarzschild_photon_sphere_radius(a::Real, x::Real) = 3.0

function kerr_equatorial_photon_sphere_radius(a::Real, x::Real)
    @assert isapprox(abs(x), 1.0; atol=1e-12) "Equatorial formula requires x = ±1"
    if x > 0
        return 2 * (1 + cos((2/3) * acos(-a)))
    else
        # x == -1
        return 2 * (1 + cos((2/3) * acos(a)))
    end
end

function kerr_polar_photon_sphere_radius(a::Real, x::Real)
    arg_denom = (1 - (a^2) / 3.0)
    @assert arg_denom > 0 "Polar formula domain violation (|a| too large for this closed form)."
    inside = (1 - a^2) / (arg_denom^(3/2))
    inside_clamped = clamp(inside, -1.0, 1.0)
    return 1.0 + 2.0 * sqrt(arg_denom) * cos((1/3) * acos(inside_clamped))
end

function kerr_extremal_photon_sphere_radius(a::Real, x::Real)
    if isapprox(a, 1.0; atol=1e-12)
        if x < sqrt(3) - 1
            return 1.0 + sqrt(2.0) * sqrt(1.0 - x) - x
        else
            return 1.0
        end
    elseif isapprox(a, -1.0; atol=1e-12)
        # use symmetry
        return kerr_extremal_photon_sphere_radius(1.0, -x)
    else
        error("kerr_extremal_photon_sphere_radius only for a ≈ ±1")
    end
end

function _u0sq_of_r(a::Real, r::Real)
    # avoid singularities: r != 1 for these formulas
    if isapprox(r, 1.0; atol=1e-14)
        throw(DomainError(r, "r too close to r = 1 (horizon in these units) for analytic formula"))
    end
    # Phi
    Phi = -((r^3 - 3.0*r^2 + a^2*r + a^2) / (a * (r - 1.0)))
    # Q
    Q = - ( r^3 * (r^3 - 6.0*r^2 + 9.0*r - 4.0*a^2) ) / ( a^2 * (r - 1.0)^2 )
    # discriminant inside sqrt
    D = (a^2 - Q - Phi^2)
    arg = D^2 + 4.0*a^2*Q
    if arg < 0
        # return complex (we'll mark as NaN to indicate invalid)
        return NaN
    end
    return (D + sqrt(arg)) / (2.0*a^2)
end


function _photon_equation(a::Real, x0::Real)
    f(r) = begin
        u0 = _u0sq_of_r(a, r)
        if isnan(u0) || !isfinite(u0)
            # cause root finder to avoid this r
            return NaN
        end
        return 1.0 - u0 - x0^2
    end
    return f
end

function kerr_geo_photon_sphere_radius_numeric(a::Real, x0::Real)
    @assert abs(x0) <= 1.0 "Inclination x must satisfy |x| ≤ 1"
    @assert !isapprox(a, 0.0; atol=1e-12) "Numeric routine only for a ≠ 0; use Schwarzschild formula for a ≈ 0"

    # analytic reference radii
    req = kerr_equatorial_photon_sphere_radius(a, sign(x0))   # equatorial of same sign
    rpolar = kerr_polar_photon_sphere_radius(a, 0.0)

    # function to root-find
    f = _photon_equation(a, x0)

    # bracket based on req,rpolar
    lo = min(req, rpolar)
    hi = max(req, rpolar)

    # sometimes the bracketed interval might be degenerate; ensure a small expansion
    if isapprox(lo, hi; atol=1e-12)
        lo = lo * 0.999
        hi = hi * 1.001
    end

    # attempt bracketed root find first
    try
        # require sign change in bracket: if there's no sign change, find_zero will fail; catch and fallback
        rroot = find_zero(f, (lo, hi))
        if !isreal(rroot)
            throw(ErrorException("Root is not real"))
        end
        # sanity: check returned r is finite and gives real-valued equation
        val = f(rroot)
        if !isfinite(val)
            throw(ErrorException("Non-finite value at root"))
        end
        return real(rroot)
    catch err1
        # fallback: try midpoint initial guess with Order1 method (secant-like)
        mid = (req + rpolar) / 2.0
        try
            rroot2 = find_zero(f, mid)   # single initial guess (method auto-chosen)
            if !isreal(rroot2)
                throw(ErrorException("Root is not real (midpoint fallback)"))
            end
            val2 = f(rroot2)
            if !isfinite(val2)
                throw(ErrorException("Non-finite value at root (midpoint fallback)"))
            end
            return real(rroot2)
        catch err2
            # final fallback: try a little outward bracket expansion and attempt again
            lo2 = lo * 0.9
            hi2 = hi * 1.1
            try
                rroot3 = find_zero(f, (lo2, hi2))
                if !isreal(rroot3)
                    throw(ErrorException("Root is not real (expanded bracket)"))
                end
                return real(rroot3)
            catch err3
                throw(ErrorException("Photon-sphere root-finding failed. Tried bracket [$(lo),$(hi)], midpoint $(mid) and expanded bracket. Errors:\n 1) $(err1)\n 2) $(err2)\n 3) $(err3)"))
            end
        end
    end
end

function kerr_geo_photon_sphere_radius(a::Real, x::Real)
    @assert abs(x) <= 1.0 "Inclination parameter x must satisfy |x| ≤ 1"

    # Schwarzschild case
    if isapprox(a, 0.0; atol=1e-12)
        return schwarzschild_photon_sphere_radius(a, x)
    end

    # Extremal analytic
    if isapprox(abs(a), 1.0; atol=1e-12)
        return kerr_extremal_photon_sphere_radius(a, x)
    end

    # Equatorial analytic
    if isapprox(abs(x), 1.0; atol=1e-12)
        return kerr_equatorial_photon_sphere_radius(a, x)
    end

    # Polar analytic (x ~ 0)
    if isapprox(abs(x), 0.0; atol=1e-14)
        return kerr_polar_photon_sphere_radius(a, x)
    end

    # Generic numeric case
    return kerr_geo_photon_sphere_radius_numeric(a, x)
end

# Innermost bound spherical orbits (IBSO)

function kerr_geo_ibso(a::Real, x::Real)
    # Schwarzschild
    if isapprox(a, 0.0; atol=1e-12)
        return 4.0
    end

    # Equatorial prograde (x = +1)
    if isapprox(x, 1.0; atol=1e-12)
        return 2 - a + 2 * sqrt(1 - a)
    end

    # Equatorial retrograde (x = -1)
    if isapprox(x, -1.0; atol=1e-12)
        return 2 + a + 2 * sqrt(1 + a)
    end

    # Polar orbits (x = 0)
    if isapprox(x, 0.0; atol=1e-12)
        δ = 27*a^4 - 8*a^6 + 3*sqrt(3) * sqrt(27*a^8 - 16*a^10)
        return 1 + sqrt(12 - 4*a^2 -
            (6*sqrt(6)*(a^2 - 2)) / sqrt(6 - 2*a^2 + 4*a^4/δ^(1/3) + δ^(1/3)) -
            4*a^4/δ^(1/3) - δ^(1/3)) / sqrt(6) +
            sqrt(6 - 2*a^2 + 4*a^4/δ^(1/3) + δ^(1/3)) / sqrt(6)
    end

    # Extremal case a=1, x=0
    if isapprox(a, 1.0; atol=1e-12) && isapprox(x, 0.0; atol=1e-12)
        return (3 + (54 - 6*sqrt(33))^(1/3) + (6*(9 + sqrt(33)))^(1/3)) / 3
    end

    # Generic case: numerical root of IBSO polynomial
    IBSOPoly(p, a, x) = ( (-4+p)^2 * p^6 +
        a^8 * (1 - x^2)^2 +
        2*a^2*p^5 * (-8 + 2*p + 4*x^2 - 3*p*x^2) +
        2*a^6*p^2 * (2 - 5*x^2 + 3*x^4) +
        a^4*p^3 * (-8*(1 - 3*x^2 + 2*x^4) + p*(6 - 14*x^2 + 9*x^4)) )

    # Initial guesses depend on sign of x
    if x >= 0
        guess1 = kerr_geo_ibso(a, 1.0)
        guess2 = kerr_geo_ibso(a, 0.0)
    else
        guess1 = kerr_geo_ibso(a, 0.0)
        guess2 = kerr_geo_ibso(a, -1.0)
    end

    return find_zero(p -> IBSOPoly(p, a, x), (guess1, guess2), Bisection())
end

# Separatrix

"""
    kerr_geo_separatrix(a, e, x)

Separatrix radius for bound geodesics in Kerr spacetime.

Special cases:
- Negative spin: symmetry a -> -a, x -> -x
- Schwarzschild (a=0): p_s = 6 + 2e
- Extremal Kerr prograde equatorial (a=1, x=1): p_s = 1 + e
- Extremal Kerr polar (a=1, x=0, e=0): analytic closed form
- Extremal Kerr polar (a=1, x=0, e=1): analytic closed form
- e=1: p_s = 2 * kerr_geo_ibso(a, x)
- General case: numerical solution of polynomial conditions
"""
function kerr_geo_separatrix(a::Real, e::Real, x::Real)
    # Negative spin symmetry
    if a < 0
        return kerr_geo_separatrix(-a, e, -x)
    end

    # Schwarzschild limit
    if isapprox(a, 0.0; atol=1e-12)
        return 6 + 2*e
    end

    # Extremal Kerr (a=1), equatorial prograde
    if isapprox(a, 1.0; atol=1e-12) && isapprox(x, 1.0; atol=1e-12)
        return 1 + e
    end

    # Extremal Kerr polar special cases
    if isapprox(a, 1.0; atol=1e-12) && isapprox(x, 0.0; atol=1e-12)
        if isapprox(e, 0.0; atol=1e-12)
            return 1 + sqrt(3) + sqrt(3 + 2*sqrt(3))
        elseif isapprox(e, 1.0; atol=1e-12)
            return (3 + (54 - 6*sqrt(33))^(1/3) + (6*(9 + sqrt(33)))^(1/3)) * (2/3)
        end
    end

    # For e=1, separatrix = 2 * IBSO
    if isapprox(e, 1.0; atol=1e-12)
        return 2 * kerr_geo_ibso(a, x)
    end

    # Define polynomials (kind of horrible)
    SepPoly(p, a, e, x) = (
        -4*(3+e)*p^11+p^12+a^12*(-1+e)^4*(1+e)^8*(-1+x)^4*(1+x)^4
        -4*a^10*(-3+e)*(-1+e)^3*(1+e)^7*p*(-1+x^2)^4
        -4*a^8*(-1+e)*(1+e)^5*p^3*(-1+x)^3*(1+x)^3*(7-7*x^2-e^2*(-13+x^2)+e^3*(-5+x^2)+7*e*(-1+x^2))
        +8*a^6*(-1+e)*(1+e)^3*p^5*(-1+x^2)^2*(3+e+12*x^2+4*e*x^2+e^3*(-5+2*x^2)+e^2*(1+2*x^2))
        -8*a^4*(1+e)^2*p^7*(-1+x)*(1+x)*(-3+e+15*x^2-5*e*x^2+e^3*(-5+3*x^2)+e^2*(-1+3*x^2))
        +4*a^2*p^9*(-7-7*e+e^3*(-5+4*x^2)+e^2*(-13+12*x^2))
        +2*a^8*(-1+e)^2*(1+e)^6*p^2*(-1+x^2)^3*(2*(-3+e)^2*(-1+x^2)
        +a^2*(e^2*(-3+x^2)-3*(1+x^2)+2*e*(1+x^2)))
        -2*p^10*(-2*(3+e)^2+a^2*(-3+6*x^2+e^2*(-3+2*x^2)+e*(-2+4*x^2)))
        +a^6*(1+e)^4*p^4*(-1+x^2)^2*(-16*(-1+e)^2*(-3-2*e+e^2)*(-1+x^2)
        +a^2*(15+6*x^2+9*x^4+e^2*(26+20*x^2-2*x^4)+e^4*(15-10*x^2+x^4)+4*e^3*(-5-2*x^2+x^4)
        -4*e*(5+2*x^2+3*x^4)))-4*a^4*(1+e)^2*p^6*(-1+x)*(1+x)*(-2*(11-14*e^2+3*e^4)*(-1+x^2)
        +a^2*(5-5*x^2-9*x^4+4*e^3*x^2*(-2+x^2)+e^4*(5-5*x^2+x^4)+e^2*(6-6*x^2+4*x^4)))
        +a^2*p^8*(-16*(1+e)^2*(-3+2*e+e^2)*(-1+x^2)+a^2*(15-36*x^2+30*x^4+e^4*(15-20*x^2+6*x^4)
        +4*e^3*(5-12*x^2+6*x^4)+4*e*(5-12*x^2+10*x^4)+e^2*(26-72*x^2+44*x^4)))
    )

    SepEquat(p, a, e) = a^4*(-3-2*e+e^2)^2 + p^2*(-6-2*e+p)^2 - 2*a^2*(1+e)*p*(14+2*e^2+3*p-e*p)
    SepPolar(p, a, e) = a^6*(e-1)^2*(1+e)^4 + p^5*(-6-2*e+p) +
                        a^2*p^3*(-4*(e-1)*(1+e)^2 + (3+e*(2+3*e))*p) -
                        a^4*(1+e)^2*p*(6+2*e^3 + 2*e*(-1+p) - 3*p - 3*e^2*(2+p))

    # Numerical cases
    if isapprox(x, 1.0; atol=1e-12)   # Equatorial prograde
        return find_zero(p -> SepEquat(p, a, e), (1+e, 6+2*e), Bisection())
    elseif isapprox(x, -1.0; atol=1e-12) # Equatorial retrograde
        return find_zero(p -> SepEquat(p, a, e), (6+2*e, 5+e+4*sqrt(1+e)), Bisection())
    elseif isapprox(x, 0.0; atol=1e-12)  # Polar
        return find_zero(p -> SepPolar(p, a, e), (1+sqrt(3)+sqrt(3+2*sqrt(3)), 8.0), Bisection())
    elseif 0 < x && x < 1  # Inclined prograde
        p1 = kerr_geo_separatrix(a, e, 1.0)
        p2 = kerr_geo_separatrix(a, e, 0.0)
        return find_zero(p -> SepPoly(p, a, e, x), (p1, p2), Bisection())
    elseif -1 < x && x < 0  # Inclined retrograde
        p1 = kerr_geo_separatrix(a, e, 0.0)
        return find_zero(p -> SepPoly(p, a, e, x), (p1, 12.0), Bisection())
    else
        error("Invalid inclination x=$x, must be in [-1,1]")
    end
end

# Innermost stable spherical orbit (ISSO)

function kerr_geo_isso(a::Real, x::Real)
    if isapprox(abs(x), 1.0, atol = 1e-12)
        return kerr_geo_isco(a, x)
    else
        return kerr_geo_separatrix(a, 0.0, x)
    end
end

"""
The next three functions

    kerr_bound_Q(a, p, e, x), kerr_scatter_Q(a, p, e, x), kerr_plunge_Q(a, p, e, x)

are to judge if the orbit with the set of parameters (a, p, e, x) is bound, scatter, or plunge orbit

"""

function kerr_bound_Q(a::Real, p::Real, e::Real, x::Real)
    ps = e > 0 ? kerr_geo_separatrix(a, e, x) : kerr_geo_ibso(a, x)
    return (p >= ps) && (0 <= e < 1)
end

function kerr_scatter_Q(a::Real, p::Real, e::Real, x::Real)
    return (p >= kerr_geo_separatrix(a, e, x)) && (e >= 1)
end

function kerr_plunge_Q(a::Real, p::Real, e::Real, x::Real)
    isBound   = kerr_bound_Q(a, p, e, x)
    isScatter = kerr_scatter_Q(a, p, e, x)
    return (!isBound) && (!isScatter)
end

"""
    kerr_geo_orbit_type(a, p, e, x) returns the general type of the orbit
"""

function kerr_geo_orbit_type(a::Real, p::Real, e::Real, x::Real)
    output = String[]

    iszero_tol(v) = isapprox(v, 0.0; atol=1e-12)

    if iszero_tol(e)
        # Circular orbits
        rph  = kerr_geo_photon_sphere_radius(a, x)
        IBSO = kerr_geo_ibso(a, x)
        ISSO = kerr_geo_isso(a, x)

        if (rph < p <= IBSO)
            output = ["Unbound", "Circular", "Unstable"]
        elseif isapprox(p, IBSO; atol=1e-12)
            output = ["MarginallyBound", "Circular", "Unstable"]
        elseif (IBSO < p < ISSO)
            output = ["Bound", "Circular", "Unstable"]
        elseif isapprox(p, ISSO; atol=1e-12)
            output = ["Bound", "Circular", "MarginallyStable"]
        elseif p > ISSO
            output = ["Bound", "Circular", "Stable"]
        elseif 0 < p <= rph
            output = ["Plunge"]
        else
            output = ["NotClassified"]
        end

        # Spherical orbit check
        if !iszero_tol(abs(x) - 1) && p > rph && !iszero_tol(a)
            push!(output, "Spherical")
        end

    else
        # Non-circular orbits
        if kerr_bound_Q(a, p, e, x) && p > 0
            output = ["Bound", "Eccentric"]
        elseif kerr_scatter_Q(a, p, e, x) && p > 0
            output = ["Scatter"]
            if iszero_tol(e - 1)
                push!(output, "Parabolic")
            elseif e > 1
                push!(output, "Hyperbolic")
            end
        elseif p > 0
            output = ["Plunge", "eccentric"]
        else
            output = ["NotClassified"]
        end
    end

    # Equatorial / Inclined classification
    if output[1] != "NotClassified" 
        if iszero_tol(abs(x) - 1)
            push!(output, "Equatorial")
        else
            push!(output, "Inclined")
        end
    end

    return output
end


end
