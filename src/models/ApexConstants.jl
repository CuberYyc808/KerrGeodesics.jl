# Constants of motion (E, Lz, Q) of the orbit with APEX parameters (a, p, e, x): one solver for
# every spin, eccentricity and inclination, and the double-root residual behind the separatrix.

# Constants of motion from the two radial roots r₁ = p/(1 − e) and r₂ = p/(1 + e), written in
# u = 1/r. With U = E², V = E L, W = L², L = Lz/x and Q = (1 − x²)(a²(1 − E²) + L²), the radial
# potential divided by r⁴ is
#
#     R/r⁴ = f U − 2x g V − h W − d,
#     f = 1 + a²(2 − x²)u² + 2a²x²u³ + a⁴(1 − x²)u⁴,   g = 2a u³,
#     h = u² − 2u³ + a²(1 − x²)u⁴,                     d = f − s,   s = 2u + 2a²u³,
#
# and it vanishes at u₁ and u₂. Schmidt's elimination of V and W leaves a quadratic in
# ν = 1 − E² whose coefficients are products of the brackets [X, Y] = X(u₁)Y(u₂) − X(u₂)Y(u₁).
# Every bracket is (u₂ − u₁) times ⟨X, Y⟩ = X(u₁) δY − Y(u₁) δX, where δX is the divided
# difference of the polynomial X between u₁ and u₂, formed from its coefficients; the common
# factor cancels from the quadratic, so nothing is lost as the two roots approach each other
# (e → 0, where ⟨X, Y⟩ → X Y′ − Y X′) and u₁ = 0 (e = 1) is an ordinary point. The
# discriminant is 16x²σ²(x²ε² + κφ) (a Plücker identity), the root of the sign of x is formed
# without cancellation, and L comes from the row at u₂, whose quadratic has the single
# positive root (D − b)/h(u₂) = c/(b + D). Neither a nor x is divided by, so a = 0, x = 0 and
# x = ±1 need no separate formulas.
# Coefficients of f, g, h, s and d = f − s in powers u⁰ … u⁴.
function _apex_coefficients(a, x)
    zm2 = 1 - x^2
    a2 = a^2
    fc = (one(a), zero(a), a2 * (2 - x^2), 2 * a2 * x^2, a2^2 * zm2)
    gc = (zero(a), zero(a), zero(a), 2 * a, zero(a))
    hc = (zero(a), zero(a), one(a), -2 * one(a), a2 * zm2)
    sc = (zero(a), 2 * one(a), zero(a), 2 * a2, zero(a))
    return (f = fc, g = gc, h = hc, s = sc, d = fc .- sc)
end

_apex_check(a, p, e, x) = (abs(a) <= 1 && p > 0 && e >= 0 && abs(x) <= 1) || throw(DomainError(
    (a, p, e, x), "APEX parameters need |a| ≤ 1, p > 0, e ≥ 0 and |x| ≤ 1."))

function _apex_constants(a::Real, p::Real, e::Real, x::Real)
    _apex_check(a, p, e, x)
    a, p, e, x = promote(float(a), float(p), float(e), float(x))
    zm2 = 1 - x^2
    a2 = a^2
    u1, u2 = (1 - e) / p, (1 + e) / p
    fc, gc, hc, sc, dc = _apex_coefficients(a, x)
    # value at u₁ and divided difference (X(u₂) − X(u₁))/(u₂ − u₁)
    q2 = u1 + u2
    q3 = u1^2 + u1 * u2 + u2^2
    q4 = q2 * (u1^2 + u2^2)
    at1(c) = c[1] + u1 * (c[2] + u1 * (c[3] + u1 * (c[4] + u1 * c[5])))
    dd(c) = c[2] + c[3] * q2 + c[4] * q3 + c[5] * q4
    f1, g1, h1, d1, s1 = at1(fc), at1(gc), at1(hc), at1(dc), at1(sc)
    δf, δg, δh, δd, δs = dd(fc), dd(gc), dd(hc), dd(dc), dd(sc)
    br(X1, δX, Y1, δY) = X1 * δY - Y1 * δX
    ρ, η, σ = br(f1, δf, h1, δh), br(f1, δf, g1, δg), br(g1, δg, h1, δh)
    τ, μ = br(s1, δs, h1, δh), br(s1, δs, g1, δg)
    κ, ε, φ = br(d1, δd, h1, δh), br(d1, δd, g1, δg), br(d1, δd, s1, δs)
    x2 = x^2
    A = ρ^2 + 4 * x2 * η * σ
    B = 2 * ρ * τ + 4 * x2 * σ * (η + μ)
    C = τ^2 + 4 * x2 * σ * μ
    radicand = x2 * ε^2 + κ * φ
    radicand >= 0 || throw(DomainError((a, p, e, x),
        "No timelike geodesic has radial roots p/(1 − e) and p/(1 + e) with this inclination."))
    root = 4 * x * sign(a) * abs(σ) * sqrt(radicand)
    # the root with +root, in the form without cancellation (root is ±0 for a = 0 and e = 1)
    ν = (B < 0) == (root < 0) ? (B + root) / (2 * A) : 2 * C / (B - root)
    # row at u₂: h₂ L² + 2x g₂ E L = s₂ − ν f₂
    du = u2 - u1
    f2, g2, h2, s2 = f1 + du * δf, g1 + du * δg, h1 + du * δh, s1 + du * δs
    c = s2 - ν * f2
    ν <= 1 || throw(DomainError((a, p, e, x),
        "No timelike geodesic has radial roots p/(1 − e) and p/(1 + e): E² = $(1 - ν) < 0."))
    E = sqrt(1 - ν)
    b = x * g2 * E
    D2 = b^2 + h2 * c
    D2 >= 0 || throw(DomainError((a, p, e, x),
        "No timelike geodesic has radial roots p/(1 − e) and p/(1 + e) with this inclination."))
    D = sqrt(D2)
    L = b >= 0 ? c / (b + D) : (D - b) / h2
    L >= 0 || throw(DomainError((a, p, e, x),
        "No timelike geodesic has radial roots p/(1 − e) and p/(1 + e): L² < 0."))
    return (E = E, Lz = x * L, Q = zm2 * (a2 * ν + L^2), L = L, ν = ν)
end

# The quotient R̂(u)/((u − u₁)(u − u₂)) of R/r⁴ = Σ cₖ uᵏ (the constants of the orbit) at u₂:
# c₂ + c₃(u₁ + 2u₂) + c₄(u₁² + 2u₁u₂ + 3u₂²) = c₄(u₂ − u₃)(u₂ − u₄) with the two remaining roots
# u₃, u₄ of R. It is zero when u₂ is a double root (the separatrix; at e = 0 the triple root of
# the ISSO), negative for a stable orbit (the next root u₃ lies above u₂) and positive inside
# the separatrix.
function _apex_separatrix_residual(a, p, e, x)
    c = _apex_constants(a, p, e, x)
    k = _apex_coefficients(a, x)
    U, V, W = c.E^2, c.E * c.L, c.L^2
    coefficient(i) = k.f[i] * U - 2 * x * k.g[i] * V - k.h[i] * W - k.d[i]
    u1, u2 = (1 - e) / p, (1 + e) / p
    return coefficient(3) + coefficient(4) * (u1 + 2 * u2) +
        coefficient(5) * (u1^2 + 2 * u1 * u2 + 3 * u2^2)
end
