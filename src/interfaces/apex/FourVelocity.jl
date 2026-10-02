# APEX reference API: four-velocity of a stable orbit given (a, p, e, x).

# =====================================================
# KerrGeodesics Mino parametrization
# =====================================================
"""
    kerr_geo_velocity_mino(a, p, e, x, initPhases, index)

Compute the four-velocity for a particle on a Kerr geodesic using Mino time parametrization.

- `a` : Kerr spin parameter
- `p` : semi-latus rectum
- `e` : eccentricity
- `x` : cosine of inclination angle
- `initPhases` : tuple of initial radial and polar phases (qr0, qθ0)
- `index` : "Contravariant" or "Covariant"
"""
function kerr_geo_velocity_mino(a, p, e, x, initPhases::Tuple{Float64,Float64}, index::String)
    # Constants of motion
    consts = kerr_geo_constants_of_motion(a, p, e, x)
    En, L, Q = consts["E"], consts["Lz"], consts["Q"]

    # Radial roots
    r1 = p / (1 - e)
    r2 = p / (1 + e)
    zm = _sqrt_nonnegative(1 - x^2)

    r3 = 1/(1 - En^2) - (r1 + r2)/2 +
        _sqrt_nonnegative((-(1/(1 - En^2)) + (r1 + r2)/2)^2 -
            (a^2 * Q)/(r1 * r2 * (1 - En^2)))
    r4 = (a^2 * Q) / (r1 * r2 * r3 * (1 - En^2))

    zp = _sqrt_nonnegative(a^2*(1 - En^2) + L^2 / (1 - zm^2))
    kr = ((r1 - r2)*(r3 - r4)) / ((r1 - r3)*(r2 - r4))
    kz = a^2*(1 - En^2)*zm^2 / zp^2

    # Fundamental frequencies
    Υr = π/(2 * Elliptic.K(kr)) *
        _sqrt_nonnegative((1 - En^2)*(r1 - r3)*(r2 - r4))
    Υθ = (π * zp)/(2 * Elliptic.K(kz))

    qr0, qθ0 = initPhases

    qr(λ) = λ * Υr + qr0
    qz(λ) = λ * Υθ + qθ0 + π / 2

    # Radial motion functions
    r(qr) = (r3*(r1 - r2)*Elliptic.Jacobi.sn(Elliptic.K(kr)/π*qr, kr)^2 - r2*(r1 - r3)) /
            ((r1 - r2)*Elliptic.Jacobi.sn(Elliptic.K(kr)/π*qr, kr)^2 - (r1 - r3))
    rprime(qr) = (2*(r1 - r2)*(r1 - r3)*(r2 - r3)*Elliptic.K(kr)*Elliptic.Jacobi.cn(qr*Elliptic.K(kr)/π, kr)*
                Elliptic.Jacobi.dn(qr*Elliptic.K(kr)/π, kr)*Elliptic.Jacobi.sn(qr*Elliptic.K(kr)/π, kr)) /
                (π*((-r1 + r3) + (r1 - r2)*Elliptic.Jacobi.sn(qr*Elliptic.K(kr)/π, kr)^2)^2)

    # Polar motion functions
    z(qθ) = zm * Elliptic.Jacobi.sn(Elliptic.K(kz)*2*qθ/π, kz)
    zprime(qθ) = (2*zm*Elliptic.K(kz)*Elliptic.Jacobi.cn(2*qθ*Elliptic.K(kz)/π, kz)*
                Elliptic.Jacobi.dn(2*qθ*Elliptic.K(kz)/π, kz))/π

    # Auxiliary functions
    Δ(qr) = r(qr)^2 + a^2 - 2*r(qr)
    Σ(qr, qθ) = r(qr)^2 + a^2 * z(qθ)^2
    Ω(qr) = sqrt(r(qr)^2 + a^2)

    # Define four-velocity components
    if index == "Contravariant"
        ut_contrav(λ) = 1/Σ(qr(λ), qz(λ)) * (Ω(qr(λ))^2 / Δ(qr(λ)) * (Ω(qr(λ))^2*En - a*L) - a^2*(1 - z(qz(λ))^2)*En + a*L)
        ur_contrav(λ) = rprime(qr(λ)) * Υr / Σ(qr(λ), qz(λ))
        uθ_contrav(λ) = -Υθ * zprime(qz(λ)) / (Σ(qr(λ), qz(λ)) * sqrt(1 - z(qz(λ))^2))
        uφ_contrav(λ) = 1/Σ(qr(λ), qz(λ)) * (a/Δ(qr(λ))*(Ω(qr(λ))^2*En - a*L) - a*En + L/(1 - z(qz(λ))^2))
        return [ut_contrav, ur_contrav, uθ_contrav, uφ_contrav]
    else
        # Covariant components
        ut_cov(λ) = -En
        ur_cov(λ) = rprime(qr(λ)) * Υr / Δ(qr(λ))
        uθ_cov(λ) = -Υθ * zprime(qz(λ)) / sqrt(1 - z(qz(λ))^2)
        uφ_cov(λ) = L
        return [ut_cov, ur_cov, uθ_cov, uφ_cov]
    end
    
end

# =====================================================
# Wrapper function
# =====================================================
"""
    kerr_geo_four_velocity(a, p, e, x; initPhases=(0.0,0.0), Covariant=false, Parametrization="Mino")

Four-velocity ``u^\\mu = dx^\\mu/d\\tau`` of the ``(a, p, e, x)`` orbit as four functions of Mino time ``\\lambda``,
`[u^t, u^r, u^θ, u^φ]`, or the covariant `[u_t, u_r, u_θ, u_φ]` with `Covariant=true`.
`initPhases = (qr0, qθ0)` are the radial and polar phases at ``\\lambda = 0``; `Parametrization` is
"Mino".
"""
function kerr_geo_four_velocity(a, p, e, x; initPhases=(0.0,0.0), Covariant=false, Parametrization="Mino")
    index = Covariant ? "Covariant" : "Contravariant"

    if Parametrization == "Mino"
        return kerr_geo_velocity_mino(a, p, e, x, initPhases, index)
    else
        error("kerr_geo_four_velocity: Parametrization must be \"Mino\" (got $(repr(Parametrization))).")
    end
end
