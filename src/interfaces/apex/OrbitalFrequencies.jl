# APEX reference API: radial/polar roots and the Mino, Boyer-Lindquist and proper-time
# frequencies of stable orbits.

_nonnegative_radicand(value) = max(0.0, real(value))
_sqrt_nonnegative(value) = sqrt(_nonnegative_radicand(value))

# [F(r+) - F(r-)] / (r+ - r-), or its |a| -> 1 limit F'((r+ + r-)/2) when the horizons merge.
function _horizon_divdiff(F, rp, rm)
    T = _float_type(rp, rm)
    # Plain divided difference while its cancellation error, ~eps/(rp - rm), stays below
    # ~1e-10 (Float64; `_tol` in T); it differs from F'(ρ0) only by F‴ (rp - rm)^2/24, so below
    # that the derivative is a symmetric difference with k Richardson steps: truncation
    # O(δ^(2k+2)) with δ = 1e-3 needs k + 1 ≥ (decimal digits)/6, two steps in Float64.
    rp - rm > _tol(T, 1e-6) && return (F(rp) - F(rm)) / (rp - rm)
    ρ0 = (rp + rm) / 2
    g(δ) = (F(ρ0 + δ) - F(ρ0 - δ)) / (2δ)
    δ = T(1e-3)
    k = cld(floor(Int, -log10(eps(T))), 6) - 1
    # weights of g(δ/2^j) extrapolated to δ = 0 in δ²: integer numerators over one denominator
    # (k = 2: (64g(δ/4) − 20g(δ/2) + g(δ))/45)
    w = [prod(big(4)^j // (big(4)^j - big(4)^i) for i in 0:k if i != j; init=big(1) // 1)
        for j in 0:k]
    D = lcm(denominator.(w))
    value = T(numerator(w[k + 1] * D)) * g(δ / 2^k)
    for j in k-1:-1:0
        value += T(numerator(w[j + 1] * D)) * g(δ / 2^j)
    end
    return value / T(D)
end

# -------------------------------------------------------------------
# Radial roots
# -------------------------------------------------------------------
"""
    kerr_geo_radial_roots(a, p, e, x; En, Lz, Q)

Return the four roots (r1, r2, r3, r4) of the radial potential, with r1 = p/(1 − e) and
r2 = p/(1 + e). `En`, `Lz` and `Q` default to the constants of motion of (a, p, e, x).
For e = 1 (E = 1) the potential is cubic: r1 = Inf, and r3, r4 are given in closed form.
"""
function kerr_geo_radial_roots(a::Real, p::Real, e::Real, x::Real; En = nothing, Lz = nothing,
        Q = nothing)
    if En === nothing || Lz === nothing || Q === nothing
        c = _apex_constants(a, p, e, x)
        En === nothing && (En = c.E)
        Lz === nothing && (Lz = c.Lz)
        Q === nothing && (Q = c.Q)
    end

    if e != 1
        r1 = p / (1 - e)
        r2 = p / (1 + e)
        # R(r1) = 0 divided by r1³ gives κ = (1 − E²) r1 and r1(2 − κ) as sums of terms that
        # shrink with 1/r1. The root sum 2/(1 − E²) = 2r1/κ then yields r3 + r4 without
        # subtracting r1 from 2/(1 − E²), which near e = 1 would leave no digits.
        w = 1 + a^2 / r1^2
        u1 = ((Lz - a * En)^2 + Q) / r1
        κ = (2 - (Lz^2 + Q) / r1 + 2u1 / r1 - a^2 * Q / r1^3) / w
        AplusB = (2a^2 / r1 + Lz^2 + Q - 2u1 + a^2 * Q / r1^2) / (w * κ) - r2
        AB = a^2 * Q / (κ * r2)

        r3 = (AplusB + _sqrt_nonnegative(AplusB^2 - 4.0 * AB)) / 2.0
        r4 = iszero(r3) ? zero(r3) : AB / r3
        return (r1, r2, r3, r4)
    end

    # Parabolic case (e == 1): use the complicated closed-form expressions
    rho2 = p / (1 + e)
    r1 = Inf
    r2 = rho2

    denom = a^2 * (-1 + x^2) - (-2 + rho2) * rho2
    inner_sqrt1 = rho2 * (a^2 + (-2 + rho2) * rho2) * (-a^2 * (-1 + x^2) + rho2^2)
    termA = 8.0 * a^2 * (-1 + x^2) * (-2.0 * a * x * rho2 + sqrt(_float_type(a, p, e, x)(2)) * _sqrt_nonnegative(inner_sqrt1))^2 / (rho2 * denom^2)
    big_inner = a^4 * (x^2 - x^4) + 2.0 * (-2 + rho2) * rho2^2 + a^2 * rho2 * (2.0 + x^2 * rho2) - 2.0 * sqrt(_float_type(a, p, e, x)(2)) * a * x * _sqrt_nonnegative(inner_sqrt1)
    termB = 4.0 * rho2^2 * big_inner^2 / denom^4
    sqrt_part = _sqrt_nonnegative(termA + termB)
    numerator = -a^4 * (-1 + x^2) * (-1 + x^2 + 2.0 * rho2) -
                4.0 * sqrt(_float_type(a, p, e, x)(2)) * a * x * rho2 * _sqrt_nonnegative(inner_sqrt1) +
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

Return `(zp, zm)` with ``z_m = \\sqrt{1 - x^2}`` and
``z_p^2 = a^2(1 - E^2) + L_z^2/x^2 = Q/z_m^2`` (``z_p = \\sqrt{Q}`` for
polar orbits, ``x = 0``); ``a^2(1 - E^2)(z_m/z_p)^2`` is the parameter of the polar elliptic functions.
"""
function kerr_geo_polar_roots(a::Real, p::Real, e::Real, x::Real)
    c = _apex_constants(a, p, e, x)
    # zp² = a²(1 − E²) + (Lz/x)² from Lz/x itself, finite at x = 0 and without 1 − z₋² for |x| ≪ 1
    return (sqrt(a^2 * c.ν + c.L^2), sqrt(1 - x^2))
end

# sign(x) with x = 0 taken from the side of its sign bit: x = 0 is the limit x → 0⁺ (a polar
# orbit's φ gains +π at every pass over the axis), and −x maps it to the limit from the other
# side, as the a < 0 reflection (a, x, ϒϕ) → (−a, −x, −ϒϕ) needs
_polar_side(x) = signbit(x) ? -one(x) : one(x)

function schwarzschild_geo_mino_frequencies(a::Real, p::Real, e::Real, x::Real)

    # Case 1: e ≈ 0
    if isapprox(e, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        return Dict(
            "ϒr" => _sqrt_nonnegative((p * (p - 6)) / (p - 3)),
            "ϒθ" => p / _sqrt_nonnegative(p - 3),
            "ϒϕ" => (p * _polar_side(x)) / _sqrt_nonnegative(p - 3),
            "ϒt" => _sqrt_nonnegative(p^5 / (p - 3))
        )
    end

    # 0 < e < 1. Toward e = 1 the characteristic n1 → 1 and ϒt ∝ (1 − e)^(−3/2) comes from
    # Π(n1|m)/(e² − 1): e² − 1 and 1 − n1 are formed from their factors.
    m = (4*e) / (p - 6 + 2*e)
    n1 = (2*e * (-4+p)) / ((1+e) * (-6+2*e+p))
    Π1 = _complete_pi(n1, (1 - e) * (p - 6 - 2e) / ((1 + e) * (p - 6 + 2e)), m)
    n2 = (16*e) / (12+8*e-4*e^2-8*p+p^2)
    Π2 = _complete_pi(n2, (p - 6 - 2e) * (p - 2 + 2e) / ((p - 6 + 2e) * (p - 2 - 2e)), m)
    e2m1 = (e - 1) * (e + 1)

    return Dict(
        "ϒr" => _sqrt_nonnegative(-(p * (-6 + 2*e + p)) / (3 + e^2 - p)) * π / (2 * _K(m)),
        "ϒθ" => p / _sqrt_nonnegative(p - 3 - e^2),
        "ϒϕ" => (p * _polar_side(x)) / _sqrt_nonnegative(p - 3 - e^2),
        "ϒt" => begin
            num = -(((-4+p) * p^2 * (-6+2*e+p) * _E(m)) / e2m1) +
                (p^2 * (28 + 4*e^2 - 12*p + p^2) * _K(m)) / e2m1 -
                (2 * (6 + 2*e - p) * (3 + e^2 - p) * p^2 * Π1) / ((-1+e) * (1+e)^2) +
                (4 * (-4+p) * p * (2 * (1+e) * _K(m) + (-6 - 2*e + p) * Π1)) / (1+e) +
                2 * (-4+p)^2 * ((-4+p) * _K(m) -
                ((6+2*e-p) * p * Π2) / (2+2*e-p))
            0.5 * _sqrt_nonnegative((-4*e^2 + (-2+p)^2) / (p * (-3 - e^2 + p))) * (8 + num / ( (-4+p)^2 * _K(m) ))
        end
    )
end

function schwarzschild_geo_boyerlindquist_frequencies(a::Real, p::Real, e::Real, x::Real)
    return Dict("Ωr" => _sqrt_nonnegative(p-6)/p^2,
        "Ωθ" => 1/p^(3/2),
        "Ωϕ" => _polar_side(x)/p^(3/2))
end

function kerr_geo_mino_frequency_r(a::Real, p::Real, e::Real, x::Real, EnLQ, roots)
    En, L, Q = EnLQ
    ρ1, ρ2, ρ3, ρ4 = roots

    if isapprox(a, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12)) && isapprox(e, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        return _sqrt_nonnegative(p*(p-6)/(p-3))
    end

    kr = ((ρ1 - ρ2) / (ρ1 - ρ3)) * ((ρ3 - ρ4) / (ρ2 - ρ4))
    # 1 − E² = 2/(ρ1 + ρ2 + ρ3 + ρ4), finite times ρ1 − ρ3 as e → 1
    return (π * _sqrt_nonnegative(2 * (ρ1 - ρ3) / (ρ1 + ρ2 + ρ3 + ρ4) * (ρ2 - ρ4))) /
        (2 * _K(kr))
end

function kerr_geo_mino_frequency_θ(a::Real, p::Real, e::Real, x::Real, EnLQ, roots)
    En, L, Q = EnLQ
    zp, zm = roots

    return π * zp / (2 * _K(a^2*(1-En^2)*(zm/zp)^2))
end

function kerr_geo_mino_frequency_ϕ(a, p, e, x, EnLQ, roots, zpzm)
    return kerr_geo_mino_frequency_ϕ_r(a, p, e, x, EnLQ, roots) +
            kerr_geo_mino_frequency_ϕ_θ(a, p, e, x, EnLQ, zpzm)
end

function kerr_geo_mino_frequency_ϕ_r(a, p, e, x, EnLQ, roots)
    En, L, Q = EnLQ
    ρ1, ρ2, ρ3, ρ4 = roots

    if isapprox(a^2, 1.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        ρin = 1 - sqrt(1 - a^2)
        ρout = 1 + sqrt(1 - a^2)
        kr = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρ4)/(ρ2-ρ4))
        hM = (ρ3 - 1) / (ρ2 - 1) * (ρ1 - ρ2) / (ρ1 - ρ3)
        return a * En * (2 / (ρ3 - 1) * (1 - (ρ2 - ρ3) / (ρ2 - 1) * _Pi(hM, kr) / _K(kr))
                + (2 - a*L / En) / (2 * (ρ3 - 1)^2) * ((2 - ((ρ1 - ρ3)*(ρ2 - ρ3)) / ((ρ1 - 1)*(ρ2 - 1))) +
                  ((ρ1 - ρ3)*(ρ2 - ρ4)*(ρ3 - 1)) / ((ρ1 - 1)*(ρ2 - 1)*(ρ4 - 1)) * _E(kr) / _K(kr) +
                  (ρ2 - ρ3) / (ρ2 - 1) * ((ρ1 - ρ3)/(ρ1 - 1) + (ρ2 - ρ3)/(ρ2 - 1) + (ρ4 - ρ3)/(ρ4 - 1) - 4) * _Pi(hM, kr) / _K(kr)))
    else
        # general Kerr
        ρin = 1 - sqrt(1 - a^2)
        ρout = 1 + sqrt(1 - a^2)
        kr = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρ4)/(ρ2-ρ4))
        Kr = _K(kr)
        hρ(ρ) = ((ρ1-ρ2)/(ρ1-ρ3)) * ((ρ3-ρ)/(ρ2-ρ))
        H(ρ) = (2*En*ρ - a*L)/(ρ3-ρ) * (1 - (ρ2-ρ3)/(ρ2-ρ) * _Pi(hρ(ρ), kr)/Kr)
        # a/(2 sqrt(1-a^2)) [H(ρout) - H(ρin)] as a divided difference (finite as |a| -> 1)
        return a * _horizon_divdiff(H, ρout, ρin)
    end
end

function kerr_geo_mino_frequency_ϕ_θ(a, p, e, x, EnLQ, zpzm)
    En, L, Q = EnLQ
    zp, zm = zpzm

    m = a^2*(1 - En^2)*(zm/zp)^2
    # Lz = 0: the limit of L Π(1 − x²|m)/K(m) is ±ϒθ, π per pass over the axis
    iszero(x) && return _polar_side(x) * π * zp / (2 * _K(m))
    # Π(z₋²|m) with 1 − z₋² = x² supplied: for |x| ≪ 1 the characteristic is within rounding of 1
    return L * _complete_pi(zm^2, x^2, m) / _K(m)
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
    # Π(hr|kr) with 1 − hr formed directly: hr → 1 as e → 1
    Πr = _complete_pi(hr, (ρ2 - ρ3)/(ρ1 - ρ3), kr)

    if isapprox(a^2, 1.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        hM = (ρ1 - ρ2)/(ρ1 - ρ3) * (ρ3 - 1)/(ρ2 - 1)
        return 5 * En + En * (0.5 * ((ρ3*(ρ1 + ρ2 + ρ3) - ρ1*ρ2) +
                (ρ1 + ρ2 + ρ3 + ρ4)*(ρ2 - ρ3) * Πr/_K(kr) +
                (ρ1 - ρ3)*(ρ2 - ρ4) * _E(kr)/_K(kr)) +
                2*(ρ3 + (ρ2 - ρ3) * Πr/_K(kr)) +
                (2*(4 - a*L/En))/(ρ3 - 1) * (1 - (ρ2 - ρ3)/(ρ2 - 1) * _Pi(hM, kr)/_K(kr)) +
                (2 - a*L/En)/(ρ3 - 1)^2 * (
                    (2 - ((ρ1 - ρ3)*(ρ2 - ρ3))/((ρ1 - 1)*(ρ2 - 1))) +
                    ((ρ1 - ρ3)*(ρ2 - ρ4)*(ρ3 - 1))/((ρ1 - 1)*(ρ2 - 1)*(ρ4 - 1)) * _E(kr)/_K(kr) +
                    (ρ2 - ρ3)/(ρ2 - 1) * ((ρ1 - ρ3)/(ρ1 - 1) + (ρ2 - ρ3)/(ρ2 - 1) + (ρ4 - ρ3)/(ρ4 - 1) - 4) * _Pi(hM, kr)/_K(kr)
                )
            )
    end
    
    term1 = (a^2 + 4) * En
    term2 = En * (0.5 * (ρ3*(ρ1 + ρ2 + ρ3) - ρ1*ρ2 + 
             (ρ1 + ρ2 + ρ3 + ρ4)*(ρ2 - ρ3) * Πr/_K(kr) +
             (ρ1 - ρ3)*(ρ2 - ρ4) * _E(kr)/_K(kr)) +
             2*(ρ3 + (ρ2 - ρ3) * Πr/_K(kr)) +
             2 * _horizon_divdiff(ρ -> ((4 - a*L/En)*ρ - 2a^2)/(ρ3 - ρ) *
                (1 - (ρ2 - ρ3)/(ρ2 - ρ) * _Pi(((ρ1 - ρ2)/(ρ1 - ρ3)) * ((ρ3 - ρ)/(ρ2 - ρ)), kr)/_K(kr)),
                ρout, ρin))
    return term1 + term2
end

function kerr_geo_mino_frequency_t_θ(a, p, e, x, params, zvals)
    En, L, Q = params
    zp, zm = zvals
    
    # E Q (1 − E(m)/K(m))/((1 − E²) z₋²) with m = a²(1 − E²)(z₋/z₊)², written with
    # D(m) = (K − E)/m so that nothing is divided by 1 − E² (→ 0 as e → 1); −a²E for Q = 0
    m = a^2*(1 - En^2)*(zm/zp)^2
    return En*Q*a^2/zp^2 * _D(m)/_K(m) - a^2*En
end

function kerr_geo_mino_frequencies(a, p, e, x)

    if isapprox(a, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        return schwarzschild_geo_mino_frequencies(a, p, e, x)
    end
    if a < 0
        freqs = kerr_geo_mino_frequencies(-a, p, e, -x)
        return Dict(
            "ϒr" => real(freqs["ϒr"]),
            "ϒθ" => real(freqs["ϒθ"]),
            "ϒϕ" => real(- freqs["ϒϕ"]),
            "ϒt" => real(freqs["ϒt"])
        )
    elseif isapprox(e, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12)) && isapprox(x, 1.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        U_r = _sqrt_nonnegative(p * (-2*a^2 + 6*a*sqrt(p) + (-5 + p)*p +
                ((a - sqrt(p))^2 * (a^2 - 4*a*sqrt(p) - (-4 + p)*p)) /
                abs(a^2 - 4*a*sqrt(p) - (-4 + p)*p)) /
                (2*a*sqrt(p) + (-3 + p)*p))
        U_theta = abs((p^(1/4) * _sqrt_nonnegative(3*a^2 - 4*a*sqrt(p) + p^2)) / _sqrt_nonnegative(2*a + (-3 + p)*sqrt(p)))
        U_phi = p^(5/4) / _sqrt_nonnegative(2*a + (-3 + p)*sqrt(p))
        U_t = (p^(5/4) * (a + p^(3/2))) / _sqrt_nonnegative(2*a + (-3 + p)*sqrt(p))

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
        r1, r2, r3, r4 = kerr_geo_radial_roots(a, p, e, x; En=En, Lz=L, Q=Q)
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
    if isapprox(a, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12)) && isapprox(e, 0.0; atol=_tol(_float_type(a, p, e, x), 1e-12))
        return schwarzschild_geo_boyerlindquist_frequencies(a, p, e, x)
    end

    MinoFreqs = kerr_geo_mino_frequencies(a, p, e, x)
    Γ = MinoFreqs["ϒt"]

    return Dict(
        "Ωr" => MinoFreqs["ϒr"] / Γ,
        "Ωθ" => MinoFreqs["ϒθ"] / Γ,
        "Ωϕ" => MinoFreqs["ϒϕ"] / Γ
    )
end

# Γτ = ⟨Σ⟩ = ⟨r²⟩ + a²⟨z²⟩ averaged over Mino time; proper-time frequency = ϒ / Γτ.
function kerr_geo_proper_frequency_factor(a, p, e, x)
    ρ1, ρ2, ρ3, ρ4 = kerr_geo_radial_roots(a, p, e, x)
    zp, zm = kerr_geo_polar_roots(a, p, e, x)
    En = kerr_geo_energy(a, p, e, x)
    kr = (ρ1 - ρ2) / (ρ1 - ρ3) * (ρ3 - ρ4) / (ρ2 - ρ4)
    kθ = a^2 * (1 - En^2) * (zm / zp)^2
    hr = (ρ1 - ρ2) / (ρ1 - ρ3)
    Kr = _K(kr)
    r2 = 0.5 * (ρ3 * (ρ1 + ρ2 + ρ3) - ρ1 * ρ2 +
                (ρ1 + ρ2 + ρ3 + ρ4) * (ρ2 - ρ3) * _complete_pi(hr, (ρ2 - ρ3) / (ρ1 - ρ3), kr) / Kr +
                (ρ1 - ρ3) * (ρ2 - ρ4) * _E(kr) / Kr)
    # a²⟨z²⟩ = a² z₋² D(kθ)/K(kθ), D = (K − E)/kθ
    a2z2 = a^2 * zm^2 * _D(kθ) / _K(kθ)
    return r2 + a2z2
end

function kerr_geo_proper_frequencies(a, p, e, x)
    f = kerr_geo_mino_frequencies(a, p, e, x)
    Γ = kerr_geo_proper_frequency_factor(a, p, e, x)
    return Dict("Ωr" => f["ϒr"] / Γ, "Ωθ" => f["ϒθ"] / Γ, "Ωϕ" => f["ϒϕ"] / Γ)
end

"""
    kerr_geo_frequencies(a, p, e, x; Time="Mino", precision=nothing)

The fundamental frequencies of the bound orbit `(a, p, e, x)`, 0 ≤ e < 1, as a `Dict`
(other eccentricities have no periodic radial motion and raise a `DomainError`). `Time="Mino"`
gives the Mino-time frequencies `"ϒr"`, `"ϒθ"`, `"ϒϕ"` and `"ϒt"` (the mean of ``dt/d\\lambda``);
`Time="BoyerLindquist"` gives `"Ωr"`, `"Ωθ"`, `"Ωϕ"`, the frequencies in coordinate time,
``\\Omega_i = \\Upsilon_i/\\Upsilon_t``, with ``i \\in \\{r, \\theta, \\phi\\}``;
`Time="Proper"` gives the frequencies in proper time, ``\\Upsilon_i`` divided by the mean
of ``d\\tau/d\\lambda``. At `x = 0` (``L_z = 0``) `"ϒϕ"` is the limit ``x \\to 0^+`` and includes
`"ϒθ"` (`x = -0.0` gives ``x \\to 0^-``). The frequencies are computed in the floating-point type
of `(a, p, e, x)`; `precision = p` converts them to `BigFloat` of `p` bits.
"""
function kerr_geo_frequencies(a, p, e, x; Time="Mino", precision=nothing)
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geo_frequencies(BigFloat(a), BigFloat(p), BigFloat(e), BigFloat(x); Time=Time)
    end
    T = _float_type(a, p, e, x)
    return _with_precision(T, _input_precision(a, p, e, x)) do
        _kerr_geo_frequencies(T(a), T(p), T(e), T(x), Time)
    end
end

function _kerr_geo_frequencies(a, p, e, x, Time)
    0 <= e < 1 || throw(DomainError(e,
        "Orbital frequencies are defined for bound orbits, 0 ≤ e < 1."))
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
