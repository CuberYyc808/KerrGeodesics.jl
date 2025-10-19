module OrbitalFrequencies

include("ConstantsOfMotion.jl")
using .ConstantsOfMotion
using Elliptic

export kerr_geo_frequencies, kerr_geo_radial_roots, kerr_geo_polar_roots

# -------------------------------------------------------------------
# Radial roots
# -------------------------------------------------------------------
"""
    kerr_geo_radial_roots(a, p, e, x; En=nothing, Q=nothing)

Return the four radial roots (r1, r2, r3, r4) of the radial equation.
If `En` or `Q` is not provided, they are computed from ConstantsOfMotion.

Two branches:
- generic bound/scattering (e != 1)
- parabolic (e == 1) uses the specialized closed-form expressions.

Notes:
- In the parabolic case r1 = Inf, r2 = ρ2 = p/(1+e).
"""
function kerr_geo_radial_roots(a::Real, p::Real, e::Real, x::Real; En = -10.0, Q = -10.0)
    # compute En and Q if not supplied
    if En == -10.0
        En = kerr_geo_energy(a, p, e, x)
    end
    if Q == -10.0
        Q = kerr_geo_carter_constant(a, p, e, x)
    end

    # Generic case (e != 1)
    if !isapprox(e, 1.0; atol=1e-12)
        r1 = p / (1 - e)
        r2 = p / (1 + e)
        AplusB = 2.0 / (1 - En^2) - (r1 + r2)  
        AB = (a^2 * Q) / ((1 - En^2) * r1 * r2) 

        r3 = (AplusB + sqrt(AplusB^2 - 4.0 * AB)) / 2.0
        r4 = AB / r3
        return (r1, r2, r3, r4)
    end

    # Parabolic case (e == 1): use the complicated closed-form expressions
    rho2 = p / (1 + e)
    r1 = Inf
    r2 = rho2

    denom = a^2 * (-1 + x^2) - (-2 + rho2) * rho2
    inner_sqrt1 = rho2 * (a^2 + (-2 + rho2) * rho2) * (-a^2 * (-1 + x^2) + rho2^2)
    termA = 8.0 * a^2 * (-1 + x^2) * (-2.0 * a * x * rho2 + sqrt(2.0) * sqrt(inner_sqrt1))^2 / (rho2 * denom^2)
    big_inner = a^4 * (x^2 - x^4) + 2.0 * (-2 + rho2) * rho2^2 + a^2 * rho2 * (2.0 + x^2 * rho2) - 2.0 * sqrt(2.0) * a * x * sqrt(inner_sqrt1)
    termB = 4.0 * rho2^2 * big_inner^2 / denom^4
    sqrt_part = sqrt(termA + termB)
    numerator = -a^4 * (-1 + x^2) * (-1 + x^2 + 2.0 * rho2) -
                4.0 * sqrt(2.0) * a * x * rho2 * sqrt(inner_sqrt1) +
                rho2^2 * (-4.0 + 4.0 * rho2 - 5.0 * rho2^2 + 2.0 * rho2^3) -
                2.0 * a^2 * rho2 * (-2.0 + 3.0 * rho2 - 2.0 * rho2^2 + x^2 * (2.0 - 5.0 * rho2 + rho2^2))
    big_frac = numerator / denom^2
    r3 = 0.25 * (1.0 - 2.0 * rho2 + sqrt_part + big_frac)
    r4 = 0.25 * (1.0 - 2.0 * rho2 - sqrt_part + big_frac)

    return (r1, r2, r3, r4)
end

# -------------------------------------------------------------------
# Polar roots
# -------------------------------------------------------------------
"""
    kerr_geo_polar_roots(a, p, e, x)

Return (zp, zm) where zm = sqrt(1 - x^2) and zp is the relevant polar amplitude.

- For generic x (not near zero), zp = sqrt[a^2 (1 - E^2) + L^2/(1 - zm^2)] where E,L come from ConstantsOfMotion.
- If x is (numerically) zero (polar orbit), use zp = sqrt(Q) instead.
"""
function kerr_geo_polar_roots(a::Real, p::Real, e::Real, x::Real)
    # Obtain constants of motion (returns a NamedTuple with fields .E, .Lz, .Q)
    consts = kerr_geo_constants_of_motion(a, p, e, x)
    En = consts["E"]
    L = consts["Lz"]
    Q = consts["Q"]

    zm = sqrt(max(0.0, 1.0 - x^2))

    if isapprox(x, 0.0; atol=1e-12)
        # polar special-case: use Q directly
        zp = sqrt(max(0.0, Q))
    else
        # generic polar amplitude
        zp = sqrt(a^2 * (1.0 - En^2) + L^2 / (1.0 - zm^2))
    end

    return (zp, zm)
end

function schwarzschild_geo_mino_frequencies(a::Real, p::Real, e::Real, x::Real)

    # Case 1: e ≈ 0
    if isapprox(e, 0.0; atol=1e-12)
        return Dict(
            "ϒr" => sqrt((p * (p - 6)) / (p - 3)),
            "ϒθ" => p / sqrt(p - 3),
            "ϒϕ" => (p * sign(x)) / sqrt(p - 3),
            "ϒt" => sqrt(p^5 / (p - 3))
        )
    end

    # Case 2: e == 1
    if isapprox(e, 1.0; atol=1e-12)
        m = (4*e) / (p - 6 + 2*e)
        return Dict(
            "ϒr" => sqrt(-(p * (-6 + 2*e + p)) / (3 + e^2 - p)) * π / (2 * Elliptic.K(m)),
            "ϒθ" => p / sqrt(p - 3 - e^2),
            "ϒϕ" => (p * sign(x)) / sqrt(p - 3 - e^2),
            "ϒt" => Inf
        )
    end

    # Case 3: 0 < e <1
    m = (4*e) / (p - 6 + 2*e)

    return Dict(
        "ϒr" => sqrt(-(p * (-6 + 2*e + p)) / (3 + e^2 - p)) * π / (2 * Elliptic.K(m)),
        "ϒθ" => p / sqrt(p - 3 - e^2),
        "ϒϕ" => (p * sign(x)) / sqrt(p - 3 - e^2),
        "ϒt" => begin
            num = -(((-4+p) * p^2 * (-6+2*e+p) * Elliptic.E(m)) / (-1+e^2)) +
                (p^2 * (28 + 4*e^2 - 12*p + p^2) * Elliptic.K(m)) / (-1+e^2) -
                (2 * (6 + 2*e - p) * (3 + e^2 - p) * p^2 * Elliptic.Π((2*e * (-4+p)) / ((1+e) * (-6+2*e+p)), π/2, m)) / ((-1+e) * (1+e)^2) +
                (4 * (-4+p) * p * (2 * (1+e) * Elliptic.K(m) + (-6 - 2*e + p) * Elliptic.Π((2*e * (-4+p)) / ((1+e) * (-6+2*e+p)), π/2, m))) / (1+e) +
                2 * (-4+p)^2 * ((-4+p) * Elliptic.K(m) -
                ((6+2*e-p) * p * Elliptic.Π((16*e) / (12+8*e-4*e^2-8*p+p^2), π/2, m)) / (2+2*e-p))
            denom = (p * (-3 - e^2 + p) * ( -4 + p)^2 )
            0.5 * sqrt((-4*e^2 + (-2+p)^2) / (p * (-3 - e^2 + p))) * (8 + num / ( (-4+p)^2 * Elliptic.K(m) ))
        end
    )
end

function schwarzschild_geo_boyerlindquist_frequencies(a::Real, p::Real, e::Real, x::Real)
    return Dict("Ωr" => sqrt(p-6)/p^2,
        "Ωθ" => 1/p^(3/2),
        "Ωϕ" => sign(x)/p^(3/2))
end

function schwarzschild_geo_proper_frequency_factor(a::Real, p::Real, e::Real, x::Real)
    if isapprox(e, 1.0; atol=1e-12)
        return Inf
    end
    if isapprox(e, 0.0; atol=1e-12)
        return p^2
    end
    m = (4*e) / (-6 + 2*e + p)
    num = (1+e) * (28 + 4*e^2 + (-12+p)*p) -
          ((1+e) * (-4+p) * (-6+2*e+p) * Elliptic.E(m) +
           2 * (6+2*e-p) * (3+e^2-p) * Elliptic.Π((2*e*(-4+p))/((1+e)*(-6+2*e+p)), π/2, m)) / Elliptic.K(m)
    denom = 2 * (-1+e) * (1+e)^2 * (-4+p)^2
    return p^2 * num / denom
end

function kerr_geo_mino_frequency_r(a::Real, p::Real, e::Real, x::Real, EnLQ, roots)
    En, L, Q = EnLQ
    ρ1, ρ2, ρ3, ρ4 = roots

    if isapprox(a, 0.0; atol=1e-12) && isapprox(e, 0.0; atol=1e-12)
        return sqrt(p*(p-6)/(p-3))
    end

    if isapprox(e, 1.0; atol=1e-12)   # e == 1
        kr = (ρ3 - ρ4) / (ρ2 - ρ4)
        return (π * sqrt(2 * (ρ2 - ρ4))) / (2 * Elliptic.K(kr))
    else
        kr = ((ρ1 - ρ2) / (ρ1 - ρ3)) * ((ρ3 - ρ4) / (ρ2 - ρ4))
        return (π * sqrt((1 - En^2) * (ρ1 - ρ3) * (ρ2 - ρ4))) / (2 * Elliptic.K(kr))
    end
end

function kerr_geo_mino_frequency_θ(a::Real, p::Real, e::Real, x::Real, EnLQ, roots)
    En, L, Q = EnLQ
    zp, zm = roots

    if isapprox(e, 1.0; atol=1e-12)   # e == 1
        return zp
    else
        return π * zp / (2 * Elliptic.K(a^2*(1-En^2)*(zm/zp)^2))
    end
end

function kerr_geo_mino_frequency_ϕ(a, p, e, x, EnLQ, roots, zpzm)
    return kerr_geo_mino_frequency_ϕ_r(a, p, e, x, EnLQ, roots) +
            kerr_geo_mino_frequency_ϕ_θ(a, p, e, x, EnLQ, zpzm)
end

function kerr_geo_mino_frequency_ϕ_r(a, p, e, x, EnLQ, roots)
    En, L, Q = EnLQ
    ρ1, ρ2, ρ3, ρ4 = roots

    if isapprox(a^2, 1.0; atol=1e-12) && isapprox(e, 1.0; atol=1e-12)
        kr = (ρ3 - ρ4) / (ρ2 - ρ4)
        hM = (ρ3 - 1) / (ρ2 - 1)
        return a * (2/(ρ3-1) * (1 - (ρ2-ρ3)/(ρ2-1) * Elliptic.Π(hM, π/2, kr)/Elliptic.K(kr)) +
                    (2 - a*L)/(2*(ρ3-1)^2) * ((2 - (ρ2-ρ3)/(ρ2-1)) + ((ρ2-ρ4)*(ρ3-1))/((ρ2-1)*(ρ4-1)) * Elliptic.E(kr)/Elliptic.K(kr) +
                    (ρ2-ρ3)/(ρ2-1) * (1 + (ρ2-ρ3)/(ρ2-1) + (ρ4-ρ3)/(ρ4-1) - 4) * Elliptic.Π(hM, π/2, kr)/Elliptic.K(kr)))
    elseif isapprox(e, 1.0; atol=1e-12)
        ρin = 1 - sqrt(1 - a^2)
        ρout = 1 + sqrt(1 - a^2)
        kr = (ρ3 - ρ4) / (ρ2 - ρ4)
        hout = (ρ3 - ρout) / (ρ2 - ρout)
        hin = (ρ3 - ρin) / (ρ2 - ρin)
        return a / (2*sqrt(1-a^2)) * ((2*ρout - a*L)/(ρ3-ρout) * (1 - (ρ2-ρ3)/(ρ2-ρout) * Elliptic.Π(hout, π/2, kr)/Elliptic.K(kr)) -
                                      (2*ρin - a*L)/(ρ3-ρin) * (1 - (ρ2-ρ3)/(ρ2-ρin) * Elliptic.Π(hin, π/2, kr)/Elliptic.K(kr)))
    elseif isapprox(a^2, 1.0; atol=1e-12)
        ρin = 1 - sqrt(1 - a^2)
        ρout = 1 + sqrt(1 - a^2)
        kr = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρ4)/(ρ2-ρ4))
        hM = (ρ3 - 1) / (ρ2 - 1) * (ρ1 - ρ2) / (ρ1 - ρ3)
        return a * En * (2 / (ρ3 - 1) * (1 - (ρ2 - ρ3) / (ρ2 - 1) * Elliptic.Π(hM, π/2, kr) / Elliptic.K(kr))
                + (2 - a*L / En) / (2 * (ρ3 - 1)^2) + ((2 - ((ρ1 - ρ3)*(ρ2 - ρ3)) / ((ρ1 - 1)*(ρ2 - 1))) +
                  ((ρ1 - ρ3)*(ρ2 - ρ4)*(ρ3 - 1)) / ((ρ1 - 1)*(ρ2 - 1)*(ρ4 - 1)) * Elliptic.E(kr) / Elliptic.K(kr) +
                  (ρ2 - ρ3) / (ρ2 - 1) * ((ρ1 - ρ3)/(ρ1 - 1) + (ρ2 - ρ3)/(ρ2 - 1) + (ρ4 - ρ3)/(ρ4 - 1) - 4) * Elliptic.Π(hM, π/2, kr) / Elliptic.K(kr)))
    else
        # general Kerr
        ρin = 1 - sqrt(1 - a^2)
        ρout = 1 + sqrt(1 - a^2)
        kr = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρ4)/(ρ2-ρ4))
        hout = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρout)/(ρ2-ρout))
        hin  = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρin)/(ρ2-ρin))
        return a/(2*sqrt(1-a^2)) * ((2*En*ρout - a*L)/(ρ3-ρout) * (1 - (ρ2-ρ3)/(ρ2-ρout) * Elliptic.Π(hout, π/2, kr)/Elliptic.K(kr)) -
                                    (2*En*ρin  - a*L)/(ρ3-ρin)  * (1 - (ρ2-ρ3)/(ρ2-ρin)  * Elliptic.Π(hin, π/2, kr)/Elliptic.K(kr)))
    end
end

function kerr_geo_mino_frequency_ϕ_θ(a, p, e, x, EnLQ, zpzm)
    En, L, Q = EnLQ
    zp, zm = zpzm

    if isapprox(x, 0.0; atol=1e-12)  # x == 0
        roots = kerr_geo_radial_roots(a, p, e, x)
        ρ1, ρ2, ρ3, ρ4 = roots
        num = (ρ1*ρ2*(a^4 + ρ1^2*ρ2^2 + a^2*((-2+ρ1)*ρ1 + (-2+ρ2)*ρ2))/2)
        den = (a^4*(2+ρ1+ρ2) + ρ1*ρ2*(ρ1^2*(-2+ρ2)+ρ1*(-2+ρ2)*ρ2-2*ρ2^2) + a^2*(ρ1^3 + ρ1^2*ρ2 + ρ1*(-4+ρ2)*ρ2 + ρ2^3))
        m = (a^2 * (a^4 + ρ1*(ρ1*(-2+ρ2)-2*ρ2)*ρ2 + a^2*(ρ1^2 + ρ2^2))) / (ρ1*ρ2*(a^4 + ρ1^2*ρ2^2 + a^2*((-2+ρ1)*ρ1 + (-2+ρ2)*ρ2)))
        return π * sqrt(num / den) / Elliptic.K(m)
    elseif isapprox(e, 1.0; atol=1e-12)
        roots = kerr_geo_radial_roots(a, p, e, x)
        ρ1, ρ2, ρ3, ρ4 = roots
        return sqrt(2) * sqrt((ρ2*(a^2 + ρ2^2)) / (a^2 + (-2 + ρ2)*ρ2))
    else
        m = a^2*(1 - En^2)*(zm/zp)^2
        return L * Elliptic.Pi(zm^2, π/2, m) / Elliptic.K(m)
    end
end

function kerr_geo_mino_frequency_t(a, p, e, x, params, rhos, zvals)
    En, L, Q = params
    ρ1, ρ2, ρ3, ρ4 = rhos
    zp, zm = zvals
    return kerr_geo_mino_frequency_t_r(a, p, e, x, params, rhos) + 
            kerr_geo_mino_frequency_t_θ(a, p, e, x, params, zvals)
end

function kerr_geo_mino_frequency_t_r(a, p, e, x, params, rhos)
    En, L, Q = params
    ρ1, ρ2, ρ3, ρ4 = rhos
    
    ρin = 1 - sqrt(1 - a^2)
    ρout = 1 + sqrt(1 - a^2)
    
    kr = (ρ1 - ρ2)/(ρ1 - ρ3) * (ρ3 - ρ4)/(ρ2 - ρ4)
    hout = (ρ1 - ρ2)/(ρ1 - ρ3) * (ρ3 - ρout)/(ρ2 - ρout)
    hin  = (ρ1 - ρ2)/(ρ1 - ρ3) * (ρ3 - ρin)/(ρ2 - ρin)
    hr   = (ρ1 - ρ2)/(ρ1 - ρ3)

    if isapprox(e, 1.0; atol=1e-12)
        return Inf
    end

    if isapprox(a^2, 1.0; atol=1e-12)
        return 5 * En + En * (0.5 * ((ρ3*(ρ1 + ρ2 + ρ3) - ρ1*ρ2) +
                (ρ1 + ρ2 + ρ3 + ρ4)*(ρ2 - ρ3) * Elliptic.Π(hr, π/2, kr)/Elliptic.K(kr) +
                (ρ1 - ρ3)*(ρ2 - ρ4) * Elliptic.E(kr)/Elliptic.K(kr)) +
                2*(ρ3 + (ρ2 - ρ3) * Elliptic.Π(hr, π/2, kr)/Elliptic.K(kr)) +
                (2*(4 - a*L/En))/(ρ3 - 1) * (1 - (ρ2 - ρ3)/(ρ2 - 1) * Elliptic.Π(hM, π/2, kr)/Elliptic.K(kr)) +
                (2 - a*L/En)/(ρ3 - 1)^2 * (
                    (2 - ((ρ1 - ρ3)*(ρ2 - ρ3))/((ρ1 - 1)*(ρ2 - 1))) +
                    ((ρ1 - ρ3)*(ρ2 - ρ4)*(ρ3 - 1))/((ρ1 - 1)*(ρ2 - 1)*(ρ4 - 1)) * Elliptic.E(kr)/Elliptic.K(kr) +
                    (ρ2 - ρ3)/(ρ2 - 1) * ((ρ1 - ρ3)/(ρ1 - 1) + (ρ2 - ρ3)/(ρ2 - 1) + (ρ4 - ρ3)/(ρ4 - 1) - 4) * Elliptic.Π(hM, π/2, kr)/Elliptic.K(kr)
                )
            )
    end
    
    term1 = (a^2 + 4) * En
    term2 = En * (0.5 * (ρ3*(ρ1 + ρ2 + ρ3) - ρ1*ρ2 + 
             (ρ1 + ρ2 + ρ3 + ρ4)*(ρ2 - ρ3) * Elliptic.Π(hr, π/2, kr)/Elliptic.K(kr) +
             (ρ1 - ρ3)*(ρ2 - ρ4) * Elliptic.E(kr)/Elliptic.K(kr)) +
             2*(ρ3 + (ρ2 - ρ3) * Elliptic.Π(hr, π/2, kr)/Elliptic.K(kr)) +
             1/sqrt(1 - a^2) * (
                ((4 - a*L/En)*ρout - 2a^2)/(ρ3 - ρout) * (1 - (ρ2 - ρ3)/(ρ2 - ρout) * Elliptic.Π(hout, π/2, kr)/Elliptic.K(kr)) -
                ((4 - a*L/En)*ρin  - 2a^2)/(ρ3 - ρin)  * (1 - (ρ2 - ρ3)/(ρ2 - ρin)  * Elliptic.Π(hin, π/2, kr)/Elliptic.K(kr))))
    return term1 + term2
end

function kerr_geo_mino_frequency_t_θ(a, p, e, x, params, zvals)
    En, L, Q = params
    zp, zm = zvals
    
    if isapprox(e, 1.0; atol=1e-12)
        return -a^2 + (a^2 * Q)/(2 * zp^2)
    elseif isapprox(x^2, 1.0; atol=1e-12)
        return -a^2 * En
    else
        return (En*Q)/((1 - En^2)*zm^2) * (1 - Elliptic.E(a^2*(1 - En^2)*(zm/zp)^2)/Elliptic.K(a^2*(1 - En^2)*(zm/zp)^2)) - a^2*En
    end
end

function kerr_geo_mino_frequencies(a, p, e, x)

    if isapprox(a, 0.0; atol=1e-12)
        return schwarzschild_geo_mino_frequencies(a, p, e, x)
    end
    if a < 0
        freqs = KerrGeoMinoFrequencies(-a, p, e, -x)
        return Dict(
            "ϒr" => real(freqs["ϒr"]),
            "ϒθ" => real(freqs["ϒθ"]),
            "ϒϕ" => real(- freqs["ϒϕ"]),
            "ϒt" => real(freqs["ϒt"])
        )
    elseif isapprox(e, 0.0; atol=1e-12) && isapprox(x, 1.0; atol=1e-12)
        U_r = sqrt(p * (-2*a^2 + 6*a*sqrt(p) + (-5 + p)*p +
                ((a - sqrt(p))^2 * (a^2 - 4*a*sqrt(p) - (-4 + p)*p)) /
                abs(a^2 - 4*a*sqrt(p) - (-4 + p)*p)) /
                (2*a*sqrt(p) + (-3 + p)*p))
        U_theta = abs((p^(1/4) * sqrt(3*a^2 - 4*a*sqrt(p) + p^2)) / sqrt(2*a + (-3 + p)*sqrt(p)))
        U_phi = p^(5/4) / sqrt(2*a + (-3 + p)*sqrt(p))
        U_t = (p^(5/4) * (a + p^(3/2))) / sqrt(2*a + (-3 + p)*sqrt(p))

        return Dict(
            "ϒr" => real(U_r),
            "ϒθ" => abs(U_theta),
            "ϒϕ" => real(U_phi),
            "ϒt" => real(U_t)
        )
    else
        consts = kerr_geo_constants_of_motion(a, p, e, x)
        En = consts["E"]
        L = consts["Lz"]
        Q = consts["Q"]
        r1, r2, r3, r4 = kerr_geo_radial_roots(a, p, e, x)
        zp, zm = kerr_geo_polar_roots(a, p, e, x)

        U_r     = kerr_geo_mino_frequency_r(a, p, e, x, [En, L, Q], [r1,r2,r3,r4])
        U_theta = kerr_geo_mino_frequency_θ(a, p, e, x, [En, L, Q], [zp, zm])
        U_phi   = kerr_geo_mino_frequency_ϕ(a, p, e, x, [En, L, Q], [r1,r2,r3,r4], [zp, zm])
        U_t     = kerr_geo_mino_frequency_t(a, p, e, x, [En, L, Q], [r1,r2,r3,r4], [zp, zm])

        return Dict(
            "ϒr" => real(U_r),
            "ϒθ" => abs(U_theta),
            "ϒϕ" => real(U_phi),
            "ϒt" => real(U_t)
        )
    end
end

function kerr_geo_boyerlindquist_frequencies(a, p, e, x)
    if isapprox(a, 0.0; atol=1e-12) && isapprox(e, 0.0; atol=1e-12)
        return schwarzschild_geo_boyerlindquist_frequencies(a, p, e, x)
    end

    if e > 1
        return Dict("Ωr"=>0, "Ωθ"=>0, "Ωϕ"=>0)
    end

    MinoFreqs = kerr_geo_mino_frequencies(a, p, e, x)
    Γ = MinoFreqs["ϒt"]

    return Dict(
        "Ωr" => MinoFreqs["ϒr"] / Γ,
        "Ωθ" => MinoFreqs["ϒθ"] / Γ,
        "Ωϕ" => MinoFreqs["ϒϕ"] / Γ
    )
end

function kerr_geo_proper_frequency_factor(a, p, e, x)
    if e >= 1
        return Inf
    end

    ρ1, ρ2, ρ3, ρ4 = kerr_geo_radial_roots(a, p, e, x)
    zp, zm = kerr_geo_polar_roots(a, p, e, x)
    T = kerr_geo_energy(a, p, e, x)

    kr = (ρ1 - ρ2)/(ρ1 - ρ3) * (ρ3 - ρ4)/(ρ2 - ρ4)
    k_theta = a^2 * (1 - T^2) * (zm / zp)^2
    hr = (ρ1 - ρ2)/(ρ1 - ρ3)

    return 0.5 * (
        -2*zp^2/( -1 + T^2 ) + ρ1*(-ρ2+ρ3) + ρ3*(ρ2+ρ3) +
        (ρ1 - ρ3)*(ρ2 - ρ4)*Elliptic.E(kr)/Elliptic.K(kr) +
        zp^2 * Elliptic.E(k_theta)/((-1 + T^2)*Elliptic.K(k_theta)) +
        (ρ2 - ρ3)*(ρ1 + ρ2 + ρ3 + ρ4)*Elliptic.Π(hr, π/2, kr)/(2*Elliptic.K(kr))
    )
end

function kerr_geo_proper_frequencies(a, p, e, x)
    if isapprox(e, 1.0; atol=1e-12)
        return Dict("Ωr"=>0, "Ωθ"=>0, "Ωϕ"=>0)
    end
    
    if isapprox(a, 0.0; atol=1e-12)
        P = schwarzschild_geo_proper_frequency_factor(a, p, e, x)
    else
        P = kerr_geo_proper_frequency_factor(a, p, e, x)
    end

    MinoFreqs = kerr_geo_mino_frequencies(a, p, e, x)

    return Dict(
        "Ωr" => MinoFreqs["ϒr"]/P,
        "Ωθ" => MinoFreqs["ϒθ"]/P,
        "Ωϕ" => MinoFreqs["ϒϕ"]/P
    )
end

function kerr_geo_frequencies(a, p, e, x; Time="Mino")
    if Time == "Mino"
        freqs = kerr_geo_mino_frequencies(a, p, e, x)
        return Dict(
            "ϒr" => real(freqs["ϒr"]),
            "ϒθ" => real(freqs["ϒθ"]),
            "ϒϕ" => real(freqs["ϒϕ"]),
            "ϒt" => real(freqs["ϒt"])
        )
    elseif Time == "BoyerLindquist"
        return kerr_geo_boyerlindquist_frequencies(a, p, e, x)
    elseif Time == "Proper"
        return kerr_geo_proper_frequencies(a, p, e, x)
    else
        error("Unknown Time option: $Time. Use \"Mino\", \"BoyerLindquist\", or \"Proper\".")
    end
end


end
