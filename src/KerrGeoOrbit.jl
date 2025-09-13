module KerrGeoOrbit

using Elliptic
using EllipticFunctions
using LinearAlgebra
using QuadGK
include("ConstantsOfMotion.jl")
using .ConstantsOfMotion
include("FourVelocity.jl")
using .FourVelocity
include("OrbitalFrequencies.jl")
using .OrbitalFrequencies
include("SpecialOrbits.jl")
using .SpecialOrbits

export kerr_geo_orbit

# ---- Kerr circular orbit (Mino time) ----

function kerr_geo_orbit_circular(a::Float64, p::Float64, e::Float64=0.0, x::Float64=1.0; initPhases=(0.0, 0.0, 0.0, 0.0))
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
    r1,r2,r3,r4 = kerr_geo_radial_roots(a, p, e, x; En = En, Q = Q)  

    # Trajectory functions (broadcastable)
    t(λ) = ((a^3*sqrt(2*a + (-3+p)*sqrt(p))*p^2 + 
            a*sqrt(2*a + (-3+p)*sqrt(p))*(-2+p)*p^3 + 
            a^2*sqrt((2*a + (-3+p)*sqrt(p))*p^7) - 
            2*sqrt((2*a + (-3+p)*sqrt(p))*p^9) + 
            sqrt((2*a + (-3+p)*sqrt(p))*p^11)) * λ) / 
            ((2*a + (-3+p)*sqrt(p))*p^(3/4)*(a^2 + (-2+p)*p))
    r(λ) = p  # constant radius
    θ(λ) = π/2 # equatorial
    ϕ(λ) = (p^(5/4)/sqrt(2*a + (-3+p)*sqrt(p))) * λ

    # Four-velocity
    velocity = kerr_geo_four_velocity(a, p, e, x; initPhases=(initPhases[2], initPhases[3]), Covariant=false, Parametrization="Mino")

    # Orbit type
    type = kerr_geo_orbit_type(a, p, e, x)

    # Associate dictionary
    return assoc = Dict(
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
        "a" => a,
        "p" => p,
        "e" => e,
        "Cosθ_inc" => x,
        "Type" => type,
        "InitialPhases" => initPhases
    )
end


function kerr_geo_orbit_generic(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))

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
    r1,r2,r3,r4 = kerr_geo_radial_roots(a, p, e, x; En, Q)
    zp, zm = kerr_geo_polar_roots(a, p, e, x)

    # Jacobi elliptic modulus for radial and polar motion
    kr = ((r1-r2)/(r1-r3)) * ((r3-r4)/(r2-r4))
    kθ = a^2 * (1-En^2) * (zm/zp)^2

    # Horizon radii
    rp = 1.0 + sqrt(1.0^2 - a^2)
    rm = 1.0 - sqrt(1.0^2 - a^2)

    # Elliptic Pi parameters for radial motion
    hr = (r1-r2)/(r1-r3)
    hp = ((r1-r2)*(r3-rp))/((r1-r3)*(r2-rp))
    hm = ((r1-r2)*(r3-rm))/((r1-r3)*(r2-rm))

    # Radial JacobiSN mapping
    rq(qr) = (r3*(r1 - r2) * Elliptic.Jacobi.sn(Elliptic.K(kr)/π * qr, kr)^2 - r2*(r1-r3)) /
            ((r1-r2) * Elliptic.Jacobi.sn(Elliptic.K(kr)/π * qr, kr)^2 - (r1-r3))

    # Polar JacobiSN mapping
    zq(qθ) = zm * Elliptic.Jacobi.sn(Elliptic.K(kθ) * 2/π * (qθ + π/2), kθ)

    # Radial and polar Jacobi amplitudes
    ψr(qr) = Elliptic.Jacobi.am(Elliptic.K(kr)/π * qr, kr)
    ψθ(qθ) = Elliptic.Jacobi.am(Elliptic.K(kθ)*2/π*(qθ+π/2), kθ)

    # t and phi increments due to radial motion

    function tr(qr)
        prefac = - En / sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
        term1 = 4 * (r2 - r3) * (Elliptic.Pi(hr, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hr, kr)))
        term2 = - 4 * (r2 - r3) / (rp - rm) * ((-1 / ((-rm + r2) * (-rm + r3))) * (-2*a^2 + rm*(4 - (a*Lz)/En)) *
            (Elliptic.Pi(hm, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hm, kr))) +
            (1 / ((-rp + r2) * (-rp + r3))) * (-2*a^2 + rp*(4 - (a*Lz)/En)) *
            (Elliptic.Pi(hp, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hp, kr))))
        term3 = (r2 - r3) * (r1 + r2 + r3 + r4) * (Elliptic.Pi(hr, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hr, kr)))
        term4 = (r1 - r3) * (r2 - r4) * (Elliptic.E(kr) * qr / π - Elliptic.E(ψr(qr), kr) +
            hr * (sin(ψr(qr)) * cos(ψr(qr)) * sqrt(1 - kr * sin(ψr(qr))^2)) / (1 - hr * sin(ψr(qr))^2))
        return prefac * (term1 + term2 + term3 + term4)
    end

    function ϕr(qr)
        prefac = 2 * a * En / ((-rm + rp) * sqrt((1 - En^2) * (r1 - r3) * (r2 - r4)))
        term_rm = (-1 / ((-rm + r2) * (-rm + r3))) * (2*rm - (a*Lz)/En) * (r2 - r3) *
            (Elliptic.Pi(hm, π/2, kr) * qr / π - real.(EllipticFunctions.ellipticPI(ψr(qr), hm, kr)))
        term_rp = (1 / ((-rp + r2) * (-rp + r3))) * (2*rp - (a*Lz)/En) * (r2 - r3) *
            (Elliptic.Pi(hp, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hp, kr)))
        return prefac * (term_rm + term_rp)
    end

    # t and phi increments due to polar motion
    tθ(qθ) = En*zp/(1-En^2) * (Elliptic.E(kθ)*2*((qθ+π/2)/π) - Elliptic.E(ψθ(qθ), kθ))
    ϕθ(qθ) = -Lz/zp * (Elliptic.Pi(zm^2, π/2, kθ)*2*((qθ+π/2)/π) - real(EllipticFunctions.ellipticPI(ψθ(qθ), zm^2, kθ)))

    qt0, qr0, qθ0, qϕ0 = initPhases

    # Total trajectory functions
    t(λ) = qt0 + ϒt * λ + tr(qr0 + ϒr * λ) + tθ(qθ0 + ϒθ * λ) - tr(qr0) - tθ(qθ0)
    r(λ) = rq(qr0 + ϒr * λ)
    θ(λ) = acos(zq(qθ0 + ϒθ * λ))
    ϕ(λ) = qϕ0 + ϒϕ * λ + ϕr(qr0 + ϒr * λ) + ϕθ(qθ0 + ϒθ * λ) - ϕr(qr0) - ϕθ(qθ0)

    # Four-velocity
    velocity = kerr_geo_four_velocity(a, p, e, x; initPhases=(initPhases[2], initPhases[3]), Covariant=false, Parametrization="Mino")

    # Orbit type
    type = kerr_geo_orbit_type(a, p, e, x)

    return assoc = Dict(
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
        "FourVelocity" => velocity,
        "Type" => type,
        "InitialPhases" => initPhases
    )
end

InverseJacobiSN(z, m) = quadgk(t -> 1/(sqrt(1-t^2)*sqrt(1-m*t^2)), 0, z)[1]

function kerr_geo_orbit_scatter(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))

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
    r1,r2,r3,r4 = kerr_geo_radial_roots(a, p, e, x; En, Q)
    zp, zm = kerr_geo_polar_roots(a, p, e, x)

    # Jacobi elliptic modulus for radial and polar motion
    kr = ((r1-r2)/(r1-r3)) * ((r3-r4)/(r2-r4))
    kθ = a^2 * (1-En^2) * (zm/zp)^2

    # Horizon radii
    rp = 1.0 + sqrt(1.0^2 - a^2)
    rm = 1.0 - sqrt(1.0^2 - a^2)

    # Elliptic Pi parameters for radial motion
    hr = (r1-r2)/(r1-r3)
    hp = ((r1-r2)*(r3-rp))/((r1-r3)*(r2-rp))
    hm = ((r1-r2)*(r3-rm))/((r1-r3)*(r2-rm))

    # Radial JacobiSN mapping
    rq(qr) = (r3*(r1 - r2) * Elliptic.Jacobi.sn(Elliptic.K(kr)/π * qr, kr)^2 - r2*(r1-r3)) /
            ((r1-r2) * Elliptic.Jacobi.sn(Elliptic.K(kr)/π * qr, kr)^2 - (r1-r3))

    # Polar JacobiSN mapping
    zq(qθ) = zm * Elliptic.Jacobi.sn(Elliptic.K(kθ) * 2/π * (qθ + π/2), kθ)

    # Radial and polar Jacobi amplitudes
    ψr(qr) = EllipticFunctions.am(Elliptic.K(kr)/π * qr, kr)
    ψθ(qθ) = EllipticFunctions.am(Elliptic.K(kθ)*2/π*(qθ+π/2), kθ)

    # t and phi increments due to radial motion
    
    function tr(qr)
        prefac = - En / sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
        term1 = 4 * (r2 - r3) * (Elliptic.Pi(hr, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hr, kr)))
        term2 = - 4 * (r2 - r3) / (rp - rm) * ((-1 / ((-rm + r2) * (-rm + r3))) * (-2*a^2 + rm*(4 - (a*Lz)/En)) *
            (Elliptic.Pi(hm, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hm, kr))) +
            (1 / ((-rp + r2) * (-rp + r3))) * (-2*a^2 + rp*(4 - (a*Lz)/En)) *
            (Elliptic.Pi(hp, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hp, kr))))
        term3 = (r2 - r3) * (r1 + r2 + r3 + r4) * (Elliptic.Pi(hr, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hr, kr)))
        term4 = (r1 - r3) * (r2 - r4) * (Elliptic.E(kr) * qr / π - Elliptic.E(ψr(qr), kr) +
            hr * (sin(ψr(qr)) * cos(ψr(qr)) * sqrt(1 - kr * sin(ψr(qr))^2)) / (1 - hr * sin(ψr(qr))^2))
        return prefac * (term1 + term2 + term3 + term4)
    end

    function ϕr(qr)
        prefac = 2 * a * En / ((-rm + rp) * sqrt((1 - En^2) * (r1 - r3) * (r2 - r4)))
        term_rm = (-1 / ((-rm + r2) * (-rm + r3))) * (2*rm - (a*Lz)/En) * (r2 - r3) *
            (Elliptic.Pi(hm, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hm, kr)))
        term_rp = (1 / ((-rp + r2) * (-rp + r3))) * (2*rp - (a*Lz)/En) * (r2 - r3) *
            (Elliptic.Pi(hp, π/2, kr) * qr / π - real(EllipticFunctions.ellipticPI(ψr(qr), hp, kr)))
        return prefac * (term_rm + term_rp)
    end

    # t and phi increments due to polar motion
    tθ(qθ) = En*zp/(1-En^2) * (Elliptic.E(kθ)*2*((qθ+π/2)/π) - Elliptic.E(ψθ(qθ), kθ))
    ϕθ(qθ) = -Lz/zp * (EllipticFunctions.ellipticPI(π/2, zm^2, kθ)*2*((qθ+π/2)/π) - real(EllipticFunctions.ellipticPI(ψθ(qθ), zm^2, kθ)))

    qrS = π * InverseJacobiSN(sqrt((r3 - r1)/(r2 - r1)), kr) / Elliptic.K(kr)
    λS = qrS / ϒr

    qt0, qr0, qθ0, qϕ0 = initPhases

    if !isapprox(qr0, 0.0; atol=1e-12) 
        println("Scattering orbits are assumed to have qr0 = 0.0")
        qr0 = 0.0
    end

    ϕS = 2 * ϒϕ * λS + ϕθ(ϒθ * λS + qθ0) - ϕθ(- ϒθ * λS + qθ0)
    θin = acos(zq(- ϒθ * λS + qθ0))
    θout = acos(zq(ϒθ * λS + qθ0))

    # Total trajectory functions
    t(λ) = qt0 + ϒt * λ + tr(qr0 + ϒr * λ) + tθ(qθ0 + ϒθ * λ) - tr(qr0) - tθ(qθ0)
    r(λ) = rq(qr0 + ϒr * λ)
    θ(λ) = acos(zq(qθ0 + ϒθ * λ))
    ϕ(λ) = qϕ0 + ϒϕ * λ + ϕr(qr0 + ϒr * λ) + ϕθ(qθ0 + ϒθ * λ) - ϕr(qr0) - ϕθ(qθ0)

    # Four-velocity
    velocity = kerr_geo_four_velocity(a, p, e, x; initPhases=(initPhases[2], initPhases[3]), Covariant=false, Parametrization="Mino")

    # Orbit type
    type = kerr_geo_orbit_type(a, p, e, x)

    return assoc = Dict(
        "a" => a,
        "p" => p,
        "e" => e,
        "Cosθ_inc" => x,
        "Parametrization" => "Mino",
        "Energy" => En,
        "AngularMomentum" => Lz,
        "CarterConstant" => Q,
        "ConstantsOfMotion" => consts,
        "angles" => Dict("ψ" => real(ϕS-π), "θin" => θin, "θout" => θout),
        "RadialRoots" => [r1, r2, r3, r4],
        "RadialFrequency" => ϒr,
        "PolarFrequency" => ϒθ,
        "AzimuthalFrequency" => ϒϕ,
        "Frequencies" => Dict("ϒt" => ϒt, "ϒr" => ϒr, "ϒθ" => ϒθ, "ϒϕ" => ϒϕ),
        "Trajectory" => [t, r, θ, ϕ],
        "FourVelocity" => velocity,
        "Type" => type,
        "InitialPhases" => initPhases
    )
end

function kerr_geo_orbit(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    if isapprox(e, 0.0; atol = 1e-12) && isapprox(abs(x), 1.0; atol = 1e-12)
        return kerr_geo_orbit_circular(a, p, e, x; initPhases = initPhases)
    elseif e > 1.0
        return kerr_geo_orbit_scatter(a, p, e, x; initPhases = initPhases)
    else
        return kerr_geo_orbit_generic(a, p, e, x; initPhases = initPhases)
    end
end

end
