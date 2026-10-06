# Plunge reference API: four-velocity of an E < 1 plunge (`generic_plunge_velocity`).

function generic_plunge_velocity_real(a, E, L, Q, (r4, r3, r2, r1), (zm, zp), (λr0, λθ0))
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    _, ξθ, kθ = _plunge_polar_parameters(a, E, L, Q)
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)

    # Radial motion functions
    r(λ) = (r3 * (r2 - r4) - r2 * (r3 - r4) * _reference_sn(ξr * (λ + λr0), kr)^2) / ((r2 - r4) - (r3 - r4) * _reference_sn(ξr * (λ + λr0), kr)^2)
    rprime(λ) = - (2 * (r2 - r3) * (r2 - r4) * (r3 - r4) * ξr * _reference_cn(ξr * (λ + λr0), kr) *
                _reference_dn(ξr * (λ + λr0), kr) * _reference_sn(ξr * (λ + λr0), kr)) /
                ((r2 - r4) - (r3 - r4) * _reference_sn(ξr * (λ + λr0), kr)^2)^2

    # Polar motion functions
    z(λ) = sqrt(zm) * _reference_sn(ξθ * (λ + λθ0), kθ)
    zprime(λ) = sqrt(zm) * ξθ * _reference_cn(ξθ * (λ + λθ0), kθ) * _reference_dn(ξθ * (λ + λθ0), kθ)

    Δ(λ) = r(λ)^2 - 2 * r(λ) + a^2
    Σ(λ) = r(λ)^2 + a^2 * z(λ)^2
    P(λ) = E * (r(λ)^2 + a^2) - a * L
    Tr(λ) = (r(λ)^2 + a^2) * P(λ) / Δ(λ)
    Tθ(λ) = - a^2 * E * (1 - z(λ)^2)
    Φr(λ) = a * P(λ) / Δ(λ)
    Φθ(λ) = L / (1 - z(λ)^2)

    ut(λ) = Tr(λ) + Tθ(λ) + a * L
    ur(λ) = rprime(λ)
    uθ(λ) = zprime(λ)
    uϕ(λ) = Φr(λ) + Φθ(λ) - a * E
    return [ut, ur, uθ, uϕ]
end

function generic_plunge_velocity_complex(a, E, L, Q, (r1, r2, A, B), (zm, zp), (λr0, λθ0))
    ξr = sqrt((1 - E^2) * A * B)
    _, ξθ, kθ = _plunge_polar_parameters(a, E, L, Q)
    kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)

    # Radial motion functions
    r(λ) = (2 * A * B * (r1 + r2) + (A - B) * (A * r2 - B * r1) * _reference_sn(ξr * (λ + λr0), kr)^2 + 
            2 * A * B * (r1 - r2) * _reference_cn(ξr * (λ + λr0), kr)) / (4 * A * B + 
            (A - B)^2 * _reference_sn(ξr * (λ + λr0), kr)^2)
    rprime(λ) = - 2 * A * B * (r1 - r2) * ξr * _reference_dn(ξr * (λ + λr0), kr) * 
                _reference_sn(ξr * (λ + λr0), kr) * (4 * A * B + 2 * (A^2 - B^2) *
                _reference_cn(ξr * (λ + λr0), kr) + 2 * (A - B)^2 * _reference_cn(ξr * (λ + λr0), kr)^2 + 
                (A - B)^2 * _reference_sn(ξr * (λ + λr0), kr)^2) / (4 * A * B + (A - B)^2 * _reference_sn(ξr * (λ + λr0), kr)^2)^2

    # Polar motion functions
    z(λ) = sqrt(zm) * _reference_sn(ξθ * (λ + λθ0), kθ)
    zprime(λ) = sqrt(zm) * ξθ * _reference_cn(ξθ * (λ + λθ0), kθ) * _reference_dn(ξθ * (λ + λθ0), kθ)

    Δ(λ) = r(λ)^2 - 2 * r(λ) + a^2
    Σ(λ) = r(λ)^2 + a^2 * z(λ)^2
    P(λ) = E * (r(λ)^2 + a^2) - a * L
    Tr(λ) = (r(λ)^2 + a^2) * P(λ) / Δ(λ)
    Tθ(λ) = - a^2 * E * (1 - z(λ)^2)
    Φr(λ) = a * P(λ) / Δ(λ)
    Φθ(λ) = L / (1 - z(λ)^2)

    ut(λ) = Tr(λ) + Tθ(λ) + a * L
    ur(λ) = rprime(λ)
    uθ(λ) = zprime(λ)
    uϕ(λ) = Φr(λ) + Φθ(λ) - a * E
    return [ut, ur, uθ, uϕ]
end

function generic_plunge_velocity_real2(a, E, L, Q, (r4, r3, r2, r1), (zm, zp), (λr0, λθ0))
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    _, ξθ, kθ = _plunge_polar_parameters(a, E, L, Q)
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    k_complete = Elliptic.K(kr)

    function r(λ)
        sn2 = _reference_sn(k_complete - ξr * (λ + λr0), kr)^2
        return (r3 * (r1 - r2) * sn2 - r2 * (r1 - r3)) /
               ((r1 - r2) * sn2 - (r1 - r3))
    end

    z(λ) = sqrt(zm) * _reference_sn(ξθ * (λ + λθ0), kθ)
    zprime(λ) = sqrt(zm) * ξθ * _reference_cn(ξθ * (λ + λθ0), kθ) * _reference_dn(ξθ * (λ + λθ0), kθ)

    Δ(λ) = r(λ)^2 - 2 * r(λ) + a^2
    P(λ) = E * (r(λ)^2 + a^2) - a * L
    Tr(λ) = (r(λ)^2 + a^2) * P(λ) / Δ(λ)
    Tθ(λ) = - a^2 * E * (1 - z(λ)^2)
    Φr(λ) = a * P(λ) / Δ(λ)
    Φθ(λ) = L / (1 - z(λ)^2)
    radial_potential(λ) = (E * (r(λ)^2 + a^2) - a * L)^2 -
                          Δ(λ) * (r(λ)^2 + (a * E - L)^2 + Q)

    ut(λ) = Tr(λ) + Tθ(λ) + a * L
    ur(λ) = -sqrt(max(radial_potential(λ), 0.0))
    uθ(λ) = zprime(λ)
    uϕ(λ) = Φr(λ) + Φθ(λ) - a * E
    return [ut, ur, uθ, uϕ]
end

"""
    generic_plunge_velocity(a, E, L, Q; initPhase=(0.0, 0.0))

Return branch-specific Mino-time four-velocity callables for an E < 1
plunge trajectory. The low-level polar component is `dz/dlambda`; the
high-level `kerr_geo_plunge` wrapper also exposes `utheta`.
"""
function generic_plunge_velocity(a, E, L, Q; initPhase = (0.0, 0.0))
    rr, cf = classify_orbit(a, E, L, Q)
    pr = polar_roots(a, E, L, Q)
    if cf == "Real1"
        return generic_plunge_velocity_real(a, E, L, Q, rr, pr, initPhase)
    elseif cf == "Complex"
        return generic_plunge_velocity_complex(a, E, L, Q, rr, pr, initPhase)
    elseif cf == "Real2"
        return generic_plunge_velocity_real2(a, E, L, Q, rr, pr, initPhase)
    else 
        @info("generic_plunge_velocity: root class $cf has no E < 1 plunge region.")
    end
end
