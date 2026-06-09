module PlungeOrbit

include("OrbitClass.jl")
using .OrbitClass
using Elliptic

export generic_plunge_orbit

function elliptic_pi(h, ψ, k)
    if h > 1
        complete = real(Elliptic.F(π/2, k) - Elliptic.Pi(k/h, π/2, k) + log(ComplexF64(-1)) / (2 * sqrt((h-k)*(h-1)/h)))
    else
        complete = Elliptic.Pi(h, π/2, k)
    end
    period = div(ψ, 1.0pi)
    remainder = abs(ψ - period * 1.0pi)
    if remainder <= 0.5pi
        if h > 1
            Π = real(Elliptic.F(remainder, k) - Elliptic.Pi(k/h, remainder, k) + log(ComplexF64((sqrt((h-k)*(h-1)/h)
            *tan(remainder) + sqrt(1-k*sin(remainder)^2))/(sqrt(1-k*sin(remainder)^2) - sqrt((h-k)*(h-1)/h)*tan(remainder)))) / (2 * sqrt((h-k)*(h-1)/h)))
        else
            Π = Elliptic.Pi(h, remainder, k)
        end
        incomplete = sign(ψ) * Π
    else
        remainder = 1.0pi - remainder
        if h > 1
            Π = real(Elliptic.F(remainder, k) - Elliptic.Pi(k/h, remainder, k) + log(ComplexF64((sqrt((h-k)*(h-1)/h)
            *tan(remainder) + sqrt(1-k*sin(remainder)^2))/(sqrt(1-k*sin(remainder)^2) - sqrt((h-k)*(h-1)/h)*tan(remainder)))) / (2 * sqrt((h-k)*(h-1)/h)))
        else
            Π = Elliptic.Pi(h, remainder, k)
        end
        incomplete = sign(ψ) * (2 * complete - Π)
    end
    if abs(period) > 0.0
        incomplete += period * complete * 2
    end
    return incomplete
end

function Atan(x, y)
    if x > 0
        return atan(y/x)
    elseif x < 0 && y >= 0
        return atan(y/x) + π
    elseif x < 0 && y < 0
        return atan(y/x) - π
    elseif x == 0 && y > 0
        return π/2
    elseif x == 0 && y < 0
        return -π/2
    else
        error("Atan is undefined for (y, x) = (0, 0)")
    end
end

function Δtθ(a, En, zm, zp)
    ξθ = sqrt(a^2 * (1 - En^2) * zp)
    kθ = zm / zp
    function tθ(λ)
        am = Elliptic.Jacobi.am(ξθ * λ, kθ)
        F = Elliptic.F(am, kθ)
        E = Elliptic.E(am, kθ)
        Iz2 = zp * (F - E) / ξθ
        I = F / ξθ
        return a^2 * En * (Iz2 - I)
    end
    return tθ
end

function Δϕθ(a, En, Lz, zm, zp)
    function ϕθ(λ)
        ξθ = sqrt(a^2 * (1 - En^2) * zp)
        kθ = zm / zp
        Π = elliptic_pi(zm, Elliptic.Jacobi.am(ξθ * λ, kθ), kθ)
        return Lz * Π / ξθ
    end
    return ϕθ
end

function real1_radial_time(a, En, Lz, roots)
    r4, r3, r2, r1 = roots
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    function tr(λ)
        ξr = sqrt((1 - En^2) * (r1 - r3) * (r2 - r4)) / 2
        kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
        hr = (r3 - r4) / (r2 - r4)
        hp = (r3 - r4) / (r2 - r4) * (r2 - rp) / (r3 - rp)
        hm = (r3 - r4) / (r2 - r4) * (r2 - rm) / (r3 - rm)
        am = Elliptic.Jacobi.am(ξr * λ, kr)
        F = Elliptic.F(am, kr)
        E = Elliptic.E(am, kr)
        Π = elliptic_pi(hr, am, kr)
        Πp = elliptic_pi(hp, am, kr)
        Πm = elliptic_pi(hm, am, kr)
        I = F / ξr
        Ir = (r2 * F - (r2 - r3) * Π) / ξr
        Ir2 = ((r2 * (r2 + r3 + r4) - r3 * r4) * F + (r1 - r3) * (r2 - r4) * E - (r2 - r3) * (r1 + r2 + r3 + r4) * Π - (r1 - r3) * (r3 - r4) * cos(am) * sin(am) * sqrt(1 - kr * sin(am)^2) / (1 - hr * sin(am)^2)) / (2 * ξr)
        Irp = ((r3 - rp) * F + (r2 - r3) * Πp) / ((r2 - rp) * (r3 - rp) * ξr)
        Irm = ((r3 - rm) * F + (r2 - r3) * Πm) / ((r2 - rm) * (r3 - rm) * ξr)
        return En * Ir2 + 2 * En * Ir + ((a^2 + 4) * En - a * Lz) * I + (((4 * En - a * Lz) * rp - 2 * a^2 * En) * Irp - ((4 * En - a * Lz) * rm - 2 * a^2 * En) * Irm) / sqrt(1 - a^2)
    end
    return tr
end

function real1_radial_phi(a, En, Lz, roots)
    r4, r3, r2, r1 = roots
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    function ϕr(λ)
        ξr = sqrt((1 - En^2) * (r1 - r3) * (r2 - r4)) / 2
        kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
        hp = (r3 - r4) / (r2 - r4) * (r2 - rp) / (r3 - rp)
        hm = (r3 - r4) / (r2 - r4) * (r2 - rm) / (r3 - rm)
        am = Elliptic.Jacobi.am(ξr * λ, kr)
        F = Elliptic.F(am, kr)
        Πp = elliptic_pi(hp, am, kr)
        Πm = elliptic_pi(hm, am, kr)
        I = F / ξr
        Irp = ((r3 - rp) * F + (r2 - r3) * Πp) / ((r2 - rp) * (r3 - rp) * ξr)
        Irm = ((r3 - rm) * F + (r2 - r3) * Πm) / ((r2 - rm) * (r3 - rm) * ξr)
        return a * ((2 * En * rp - a * Lz) * Irp - (2 * En * rm - a * Lz) * Irm) / (2 * sqrt(1 - a^2)) + a * En * I 
    end
    return ϕr
end

function Δtr_complex(a, En, Lz, roots)
    r1, r2, A, B = roots
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    function tr(λ)
        ξr = sqrt((1 - En^2) * A * B)
        kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)
        f = 4 * A * B / (A - B)^2
        Dp = - sqrt(4 * A * B * (r1 - rp) * (rp - r2)) / (A * (rp - r2) + B * (r1 - rp))
        Dm = - sqrt(4 * A * B * (r1 - rm) * (rm - r2)) / (A * (rm - r2) + B * (r1 - rm))
        x = Elliptic.Jacobi.am(ξr * λ, kr)
        F = Elliptic.F(x, kr)
        E = Elliptic.E(x, kr)
        Πf = elliptic_pi(-1/f, x, kr)
        Πp = elliptic_pi(1/Dp^2, x, kr)
        Πm = elliptic_pi(1/Dm^2, x, kr)
        Ir = (A * r2 - B * r1) * λ / (A - B) + (A + B) * (r1 - r2) * Πf / (2 * (A - B) * ξr) + atan((r1 - r2) * sin(x) / sqrt(4 * A * B * (1 - kr * sin(x)^2))) / sqrt(1 - En^2) 
        Ir2 = (A * r2^2 - B * r1^2) * λ / (A - B) + sqrt(A * B / (1 - En^2)) * E - (A + B) * (A^2 + 2 * r2^2 - B^2 - 2 * r1^2) * Πf / (4 * (A - B) * ξr) + sqrt(A * B / (1 - En^2)) * ((A + B) / (A - B) + cos(x)) * sin(x) * sqrt(1 -  kr * sin(x)^2) / (f + sin(x)^2) - (A^2 + 2 * r2^2 - B^2 - 2 * r1^2) * Atan(f - (1 + 2 * f * kr) * sin(x)^2, 2 * sin(x) * sqrt(f * (1 - kr * sin(x)^2) * (1 + f * kr))) / (4 * (r1 - r2) * sqrt(1 - En^2)) 
        Irp = (A - B) * λ / (A * (r2 - rp) - B * (r1 - rp)) - (r1 - r2) * (B * (r1 - rp) - A * (rp - r2)) * Πp / (2 * ξr * (r1 - rp) * (rp - r2) * (B * (r1 - rp) + A * (rp - r2))) - sqrt((r1 - r2)/((1 - En^2) * (r1 - rp) * (rp - r2) * (A^2 * (rp - r2) + B^2 * (r1 - rp) - (r1 - r2) * (r1 - rp) * (rp - r2)))) * log(((Dp * sqrt(1 - Dp^2 * kr) + sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dp^2 - sin(x)^2)^2) / ((Dp * sqrt(1 - Dp^2 * kr) - sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dp^2 - sin(x)^2)^2)) / 4 
        Irm = (A - B) * λ / (A * (r2 - rm) - B * (r1 - rm)) - (r1 - r2) * (B * (r1 - rm) - A * (rm - r2)) * Πm / (2 * ξr * (r1 - rm) * (rm - r2) * (B * (r1 - rm) + A * (rm - r2))) - sqrt((r1 - r2)/((1 - En^2) * (r1 - rm) * (rm - r2) * (A^2 * (rm - r2) + B^2 * (r1 - rm) - (r1 - r2) * (r1 - rm) * (rm - r2)))) * log(((Dm * sqrt(1 - Dm^2 * kr) + sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dm^2 - sin(x)^2)^2) / ((Dm * sqrt(1 - Dm^2 * kr) - sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dm^2 - sin(x)^2)^2)) / 4 
        return En * Ir2 + 2 * En * Ir + ((a^2 + 4) * En - a * Lz) * λ + (((4 * En - a * Lz) * rp - 2 * a^2 * En) * Irp - ((4 * En - a * Lz) * rm - 2 * a^2 * En) * Irm) / sqrt(1 - a^2)
    end
    return tr
end

function Δϕr_complex(a, En, Lz, roots)
    r1, r2, A, B = roots
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    function ϕr(λ)
        ξr = sqrt((1 - En^2) * A * B)
        kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)
        Dp = - sqrt(4 * A * B * (r1 - rp) * (rp - r2)) / (A * (rp - r2) + B * (r1 - rp))
        Dm = - sqrt(4 * A * B * (r1 - rm) * (rm - r2)) / (A * (rm - r2) + B * (r1 - rm))
        x = Elliptic.Jacobi.am(ξr * λ, kr)
        Πp = elliptic_pi(1/Dp^2, x, kr)
        Πm = elliptic_pi(1/Dm^2, x, kr)
        Irp = (A - B) * λ / (A * (r2 - rp) - B * (r1 - rp)) - (r1 - r2) * (B * (r1 - rp) - A * (rp - r2)) * Πp / (2 * ξr * (r1 - rp) * (rp - r2) * (B * (r1 - rp) + A * (rp - r2))) - sqrt((r1 - r2)/((1 - En^2) * (r1 - rp) * (rp - r2) * (A^2 * (rp - r2) + B^2 * (r1 - rp) - (r1 - r2) * (r1 - rp) * (rp - r2)))) * log(((Dp * sqrt(1 - Dp^2 * kr) + sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dp^2 - sin(x)^2)^2) / ((Dp * sqrt(1 - Dp^2 * kr) - sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dp^2 - sin(x)^2)^2)) / 4 
        Irm = (A - B) * λ / (A * (r2 - rm) - B * (r1 - rm)) - (r1 - r2) * (B * (r1 - rm) - A * (rm - r2)) * Πm / (2 * ξr * (r1 - rm) * (rm - r2) * (B * (r1 - rm) + A * (rm - r2))) - sqrt((r1 - r2)/((1 - En^2) * (r1 - rm) * (rm - r2) * (A^2 * (rm - r2) + B^2 * (r1 - rm) - (r1 - r2) * (r1 - rm) * (rm - r2)))) * log(((Dm * sqrt(1 - Dm^2 * kr) + sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dm^2 - sin(x)^2)^2) / ((Dm * sqrt(1 - Dm^2 * kr) - sqrt(1 - kr * sin(x)^2) * sin(x))^2 + kr * (Dm^2 - sin(x)^2)^2)) / 4 
        return a * ((2 * En * rp - a * Lz) * Irp - (2 * En * rm - a * Lz) * Irm) / (2 * sqrt(1 - a^2)) + a * En * λ
    end
    return ϕr
end

function real2_radial_position(absλ, E, roots)
    r4, r3, r2, r1 = roots
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    sn2 = Elliptic.Jacobi.sn(Elliptic.K(kr) - ξr * absλ, kr)^2
    return (r3 * (r1 - r2) * sn2 - r2 * (r1 - r3)) /
           ((r1 - r2) * sn2 - (r1 - r3))
end

function real2_radial_elliptic_pi_continuation(h, ψ, k)
    if h <= 1
        return Elliptic.Pi(h, ψ, k)
    end
    α = sqrt((h - k) * (h - 1) / h)
    arg = (α * tan(ψ) + sqrt(1 - k * sin(ψ)^2)) /
          (sqrt(1 - k * sin(ψ)^2) - α * tan(ψ))
    return real(Elliptic.F(ψ, k) - Elliptic.Pi(k / h, ψ, k) + log(ComplexF64(arg)) / (2 * α))
end

function real2_radial_elliptic_pi_complete(h, k)
    h <= 1 && return Elliptic.Pi(h, π / 2, k)
    return real(Elliptic.K(k) - Elliptic.Pi(k / h, π / 2, k))
end

function real2_radial_time_phi(a, E, L, roots, λr0)
    r4, r3, r2, r1 = roots
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    K = Elliptic.K(kr)
    hr = (r1 - r2) / (r1 - r3)
    hp = ((r1 - r2) * (r3 - rp)) / ((r1 - r3) * (r2 - rp))
    hm = ((r1 - r2) * (r3 - rm)) / ((r1 - r3) * (r2 - rm))

    q_of_absλ(absλ) = π - π * ξr * absλ / K
    ψr(q) = Elliptic.Jacobi.am(K / π * q, kr)
    Πr(h, q) = real2_radial_elliptic_pi_complete(h, kr) * q / π -
               real2_radial_elliptic_pi_continuation(h, ψr(q), kr)

    function raw_tr(absλ)
        q = q_of_absλ(absλ)
        ψ = ψr(q)
        prefac = -E / sqrt((1 - E^2) * (r1 - r3) * (r2 - r4))
        term1 = 4 * (r2 - r3) * Πr(hr, q)
        term2 = -4 * (r2 - r3) / (rp - rm) *
            ((-1 / ((-rm + r2) * (-rm + r3))) * (-2 * a^2 + rm * (4 - (a * L) / E)) * Πr(hm, q) +
             (1 / ((-rp + r2) * (-rp + r3))) * (-2 * a^2 + rp * (4 - (a * L) / E)) * Πr(hp, q))
        term3 = (r2 - r3) * (r1 + r2 + r3 + r4) * Πr(hr, q)
        term4 = (r1 - r3) * (r2 - r4) *
            (Elliptic.E(kr) * q / π - Elliptic.E(ψ, kr) +
             hr * (sin(ψ) * cos(ψ) * sqrt(1 - kr * sin(ψ)^2)) / (1 - hr * sin(ψ)^2))
        return -prefac * (term1 + term2 + term3 + term4)
    end

    function raw_phir(absλ)
        q = q_of_absλ(absλ)
        prefac = 2 * a * E / ((rp - rm) * sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)))
        term_rm = (-1 / ((-rm + r2) * (-rm + r3))) * (2 * rm - (a * L) / E) * (r2 - r3) * Πr(hm, q)
        term_rp = (1 / ((-rp + r2) * (-rp + r3))) * (2 * rp - (a * L) / E) * (r2 - r3) * Πr(hp, q)
        return -prefac * (term_rm + term_rp)
    end

    dψr(q) = Elliptic.Jacobi.dn(K / π * q, kr) * K / π
    dΠr_dq(h, q) = real2_radial_elliptic_pi_complete(h, kr) / π -
                   dψr(q) / ((1 - h * sin(ψr(q))^2) * sqrt(1 - kr * sin(ψr(q))^2))

    function d_raw_tr_dq(q)
        ψ = ψr(q)
        dψ = dψr(q)
        dΠhr = dΠr_dq(hr, q)
        dΠhp = dΠr_dq(hp, q)
        dΠhm = dΠr_dq(hm, q)
        prefac = -E / sqrt((1 - E^2) * (r1 - r3) * (r2 - r4))
        d_term4 = (r1 - r3) * (r2 - r4) *
            (Elliptic.E(kr) / π -
             hr * kr * cos(ψ)^2 * sin(ψ)^2 * dψ / ((1 - hr * sin(ψ)^2) * sqrt(1 - kr * sin(ψ)^2)) -
             sqrt(1 - kr * sin(ψ)^2) * dψ +
             2 * hr^2 * cos(ψ)^2 * sin(ψ)^2 * sqrt(1 - kr * sin(ψ)^2) * dψ / (1 - hr * sin(ψ)^2)^2 +
             hr * cos(ψ)^2 * sqrt(1 - kr * sin(ψ)^2) * dψ / (1 - hr * sin(ψ)^2) -
             hr * sin(ψ)^2 * sqrt(1 - kr * sin(ψ)^2) * dψ / (1 - hr * sin(ψ)^2))
        d_terms = 4 * (r2 - r3) * dΠhr +
                  (r2 - r3) * (r1 + r2 + r3 + r4) * dΠhr -
                  4 * (r2 - r3) / (rp - rm) *
                  ((-1 / ((-rm + r2) * (-rm + r3))) * (-2 * a^2 + rm * (4 - (a * L) / E)) * dΠhm +
                   (1 / ((-rp + r2) * (-rp + r3))) * (-2 * a^2 + rp * (4 - (a * L) / E)) * dΠhp) +
                  d_term4
        return -prefac * d_terms
    end

    function d_raw_phir_dq(q)
        dΠhp = dΠr_dq(hp, q)
        dΠhm = dΠr_dq(hm, q)
        prefac = 2 * a * E / ((rp - rm) * sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)))
        d_terms = (-1 / ((-rm + r2) * (-rm + r3))) * (2 * rm - (a * L) / E) * (r2 - r3) * dΠhm +
                  (1 / ((-rp + r2) * (-rp + r3))) * (2 * rp - (a * L) / E) * (r2 - r3) * dΠhp
        return -prefac * d_terms
    end

    dq_dabsλ = -π * ξr / K
    d_raw_tr_dabsλ(absλ) = d_raw_tr_dq(q_of_absλ(absλ)) * dq_dabsλ
    d_raw_phir_dabsλ(absλ) = d_raw_phir_dq(q_of_absλ(absλ)) * dq_dabsλ

    function radial_tr_prime(absλ)
        r = real2_radial_position(absλ, E, roots)
        Δ = r^2 - 2 * r + a^2
        P = E * (r^2 + a^2) - a * L
        return ((r^2 + a^2) * P) / Δ
    end

    function radial_phir_prime(absλ)
        r = real2_radial_position(absλ, E, roots)
        Δ = r^2 - 2 * r + a^2
        P = E * (r^2 + a^2) - a * L
        return a * P / Δ
    end

    tr_linear = radial_tr_prime(λr0) - d_raw_tr_dabsλ(λr0)
    phir_linear = radial_phir_prime(λr0) - d_raw_phir_dabsλ(λr0)

    tr(absλ) = raw_tr(absλ) + tr_linear * absλ
    phir(absλ) = raw_phir(absλ) + phir_linear * absλ
    dtr(absλ) = d_raw_tr_dabsλ(absλ) + tr_linear
    dphir(absλ) = d_raw_phir_dabsλ(absλ) + phir_linear
    return tr, phir, (tr_linear=tr_linear, phir_linear=phir_linear, dtr=dtr, dphir=dphir)
end

function real2_radial_lambda_of_r(a, E, roots, r)
    r4, r3, r2, r1 = roots
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    u = (r - r2) * (r1 - r3) / ((r1 - r2) * (r - r3))
    u = clamp(u, 0.0, 1.0)
    return (Elliptic.K(kr) - Elliptic.F(asin(sqrt(u)), kr)) / ξr
end

"""
    generic_plunge_orbit(a, E, L, Q; initPhases=(0.0, 0.0, 0.0, 0.0))

Return branch-specific callable Boyer-Lindquist trajectory functions
`t(lambda)`, `r(lambda)`, `theta(lambda)`, and `phi(lambda)` for a bound plunge
Kerr plunge where the classified branch is supported.
"""
function generic_plunge_orbit(a, E, L, Q; initPhases = (0.0, 0.0, 0.0, 0.0), real2_horizon_offset=1e-4)
    roots, cf = classify_orbit(a, E, L, Q)
    zm, zp = polar_roots(a, E, L, Q)
    λt0, λr0, λθ0, λϕ0 = initPhases
    ξθ = sqrt(a^2 * (1 - E^2) * zp)
    kθ = zm / zp

    if cf == "Real1"
        r4, r3, r2, r1 = roots
        ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
        kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
        r_real1(λ) = (r3 * (r2 - r4) - r2 * (r3 - r4) * Elliptic.Jacobi.sn(ξr * (λ + λr0), kr)^2) / ((r2 - r4) - (r3 - r4) * Elliptic.Jacobi.sn(ξr * (λ + λr0), kr)^2)
        θ_real1(λ) = acos(sqrt(zm) * Elliptic.Jacobi.sn(ξθ * (λ + λθ0), kθ))
        tr_real1 = real1_radial_time(a, E, L, roots)
        ϕr_real1 = real1_radial_phi(a, E, L, roots)
        tθ_real1 = Δtθ(a, E, zm, zp)
        ϕθ_real1 = Δϕθ(a, E, L, zm, zp)
        t_real1(λ) = tr_real1(λ + λr0) - tr_real1(λr0) + tθ_real1(λ + λθ0) - tθ_real1(λθ0) + a * L * λ + λt0
        ϕ_real1(λ) = ϕr_real1(λ + λr0) - ϕr_real1(λr0) + ϕθ_real1(λ + λθ0) - ϕθ_real1(λθ0) - a * E * λ + λϕ0
        return [t_real1, r_real1, θ_real1, ϕ_real1]
    elseif cf == "Complex"
        r1, r2, A, B = roots
        ξr = sqrt((1 - E^2) * A * B)
        kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)
        r_complex(λ) = (2 * A * B * (r1 + r2) + (A - B) * (A * r2 - B * r1) * Elliptic.Jacobi.sn(ξr * (λ + λr0), kr)^2 + 
            2 * A * B * (r1 - r2) * Elliptic.Jacobi.cn(ξr * (λ + λr0), kr)) / (4 * A * B + 
            (A - B)^2 * Elliptic.Jacobi.sn(ξr * (λ + λr0), kr)^2)
        θ_complex(λ) = acos(sqrt(zm) * Elliptic.Jacobi.sn(ξθ * (λ + λθ0), kθ))
        tr_complex = Δtr_complex(a, E, L, roots)
        ϕr_complex = Δϕr_complex(a, E, L, roots)
        tθ_complex = Δtθ(a, E, zm, zp)
        ϕθ_complex = Δϕθ(a, E, L, zm, zp)
        t_complex(λ) = tr_complex(λ + λr0) - tr_complex(λr0) + tθ_complex(λ + λθ0) - tθ_complex(λθ0) + a * L * λ + λt0
        ϕ_complex(λ) = ϕr_complex(λ + λr0) - ϕr_complex(λr0) + ϕθ_complex(λ + λθ0) - ϕθ_complex(λθ0) - a * E * λ + λϕ0
        return [t_complex, r_complex, θ_complex, ϕ_complex]
    elseif cf == "Real2"
        r4, r3, r2, r1 = roots
        r_real2(λ) = real2_radial_position(λ + λr0, E, roots)
        θ_real2(λ) = acos(sqrt(zm) * Elliptic.Jacobi.sn(ξθ * (λ + λθ0), kθ))
        tr_real2, ϕr_real2, _ = real2_radial_time_phi(a, E, L, roots, λr0)
        tθ_real2 = Δtθ(a, E, zm, zp)
        ϕθ_real2 = Δϕθ(a, E, L, zm, zp)
        t_real2(λ) = tr_real2(λ + λr0) - tr_real2(λr0) + tθ_real2(λ + λθ0) - tθ_real2(λθ0) + a * L * λ + λt0
        ϕ_real2(λ) = ϕr_real2(λ + λr0) - ϕr_real2(λr0) + ϕθ_real2(λ + λθ0) - ϕθ_real2(λθ0) - a * E * λ + λϕ0
        return [t_real2, r_real2, θ_real2, ϕ_real2]
    else 
        @info("The orbit is classified as $cf which has not been included in the implementation.")
    end
end


end
