# Elliptic integrals and Jacobi functions shared by the radial models, the polar engine and the
# APEX interfaces: Carlson's R_F, R_C, R_D and R_J by duplication (Cauchy principal value for
# R_J with a negative fourth argument), the Legendre forms F, E, D = (F − E)/m, Π and the pole
# primitive (Π − F)/n built on them, the complete integrals, the Mino time to infinity in
# Carlson form, and the Landen (AGM) sequence with Bulirsch's recursion for sn, cn, dn. Every
# Legendre form takes the complement m1 = 1 − m (and n1 = 1 − n for Π) as the caller forms it
# from its roots, so a parameter within rounding of 1 keeps its digits.

# Carlson's R_F(x, y, z) by duplication (x, y, z ≥ 0, at most one of them zero, or x ≥ 0 and
# y, z complex conjugates)
function _carlson_rf(x, y, z)
    A0 = (x + y + z) / 3
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (3eps())^(1 / 6)
    xn, yn, zn, A, scale = x, y, z, A0, 1.0
    while scale * Q > abs(A)
        sx, sy, sz = sqrt(xn), sqrt(yn), sqrt(zn)
        λ = sx * sy + sx * sz + sy * sz
        xn, yn, zn, A = (xn + λ) / 4, (yn + λ) / 4, (zn + λ) / 4, (A + λ) / 4
        scale /= 4
    end
    X = (A0 - x) * scale / A
    Y = (A0 - y) * scale / A
    Z = -(X + Y)
    E2 = X * Y - Z^2
    E3 = X * Y * Z
    return (1 - E2 / 10 + E3 / 14 + E2^2 / 24 - 3 * E2 * E3 / 44) / sqrt(A)
end

# Carlson's R_C(x, y) = R_F(x, y, y); for y < 0 the Cauchy principal value
# R_C(x, y) with ∂R_C/∂y: for y > 0 by R_F(x, y, y), for y < 0 the principal value
# √(x/(x − y)) R_C(x − y, −y). Next to y = x the series R_C(x, x(1 + ε)) = x^(−1/2) Σ (−ε)^k/(2k + 1)
# gives both (the derivative's closed form (R_C − √x/y)/(2(x − y)) cancels there).
function _carlson_rc_dy(x, y)
    if y < 0
        u, v = x - y, -y
        value, dv = _carlson_rc_dy(u, v)
        du = (-value / 2 - v * dv) / u             # homogeneity: u ∂_u + v ∂_v = −R_C/2
        return sqrt(x / u) * value, sqrt(x) * (value / (2 * u^1.5) - (du + dv) / sqrt(u))
    end
    x > 0 || return pi / (2 * sqrt(y)), -pi / (4 * y^1.5)
    ε = (y - x) / x
    if abs(ε) < 0.25
        # term = (−ε)^k, dterm = (−1)^k ε^(k − 1)
        value = 0.0; derivative = 0.0; term = 1.0; dterm = -1.0
        for k in 0:32
            value += term / (2k + 1)
            k > 0 && (derivative += k * dterm / (2k + 1))
            k > 0 && abs(term) < 1e-18 && break
            term *= -ε
            k > 0 && (dterm *= -ε)
        end
        return value / sqrt(x), derivative / x^1.5
    end
    value = _carlson_rf(x, y, y)
    return value, (value - sqrt(x) / y) / (2 * (x - y))
end

# Carlson's R_D(x, y, z) by duplication (x, y ≥ 0, at most one of them zero, z > 0)
function _carlson_rd(x, y, z)
    A0 = (x + y + 3z) / 5
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (eps() / 4)^(1 / 6)
    xn, yn, zn, A, scale, total = x, y, z, A0, 1.0, 0.0
    while scale * Q > abs(A)
        sx, sy, sz = sqrt(xn), sqrt(yn), sqrt(zn)
        λ = sx * sy + sx * sz + sy * sz
        total += scale / (sz * (zn + λ))
        xn, yn, zn, A = (xn + λ) / 4, (yn + λ) / 4, (zn + λ) / 4, (A + λ) / 4
        scale /= 4
    end
    X = (A0 - x) * scale / A
    Y = (A0 - y) * scale / A
    Z = -(X + Y) / 3
    E2 = X * Y - 6 * Z^2
    E3 = (3 * X * Y - 8 * Z^2) * Z
    E4 = 3 * (X * Y - Z^2) * Z^2
    E5 = X * Y * Z^3
    return scale * (1 - 3 * E2 / 14 + E3 / 6 + 9 * E2^2 / 88 - 3 * E4 / 22 - 9 * E2 * E3 / 52 +
        3 * E5 / 26) / (A * sqrt(A)) + 3 * total
end

# Carlson's R_J(x, y, z, p) with ∂R_J/∂p, the derivative carried through the duplication itself
# (its closed form divides by (p − x)(p − y)(p − z)). For p < 0 the principal value of Carlson
# (1995) 2.27; the value alone is `_carlson_rj`.
function _carlson_rj_dp(x, y, z, p)
    if p <= 0
        xt, zt = min(x, y, z), max(x, y, z)
        yt = x + y + z - xt - zt
        a = 1 / (yt - p); da = a^2
        b = a * (zt - yt) * (yt - xt); db = da * (zt - yt) * (yt - xt)
        pt = yt + b
        rj, drj = _carlson_rj_dp(xt, yt, zt, pt)
        rc, drc = _carlson_rc_dy(xt * zt / yt, p * pt / yt)
        rf = _carlson_rf(xt, yt, zt)
        return a * (b * rj + 3 * (rc - rf)),
            da * (b * rj + 3 * (rc - rf)) + a * (db * rj + b * drj * db + 3 * drc * (pt + p * db) / yt)
    end
    A0 = (x + y + z + 2p) / 5
    δ = (p - x) * (p - y) * (p - z)
    dδ = (p - y) * (p - z) + (p - x) * (p - z) + (p - x) * (p - y)
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z), abs(A0 - p)) / (eps() / 4)^(1 / 6)
    xn, yn, zn, pn, A = x, y, z, p, A0
    dxn = dyn = dzn = 0.0; dpn = 1.0; dA = 2 / 5
    scale, total, dtotal = 1.0, 0.0, 0.0
    while scale * Q > abs(A)
        sx, sy, sz, sp = sqrt(xn), sqrt(yn), sqrt(zn), sqrt(pn)
        # an argument that is zero (a complete integral) has no p dependence at that step
        dsx, dsy, dsz = (iszero(dxn) ? 0.0 : dxn / (2sx)), (iszero(dyn) ? 0.0 : dyn / (2sy)),
            (iszero(dzn) ? 0.0 : dzn / (2sz))
        dsp = dpn / (2sp)
        λ = sx * sy + sx * sz + sy * sz
        dλ = dsx * (sy + sz) + dsy * (sx + sz) + dsz * (sx + sy)
        d = (sp + sx) * (sp + sy) * (sp + sz)
        dd = (dsp + dsx) * (sp + sy) * (sp + sz) + (sp + sx) * (dsp + dsy) * (sp + sz) +
            (sp + sx) * (sp + sy) * (dsp + dsz)
        arg = 1 + scale^3 * δ / d^2
        darg = scale^3 * (dδ / d^2 - 2 * δ * dd / d^3)
        rc, drc = _carlson_rc_dy(1.0, arg)
        total += scale * rc / d
        dtotal += scale * (drc * darg / d - rc * dd / d^2)
        xn, yn, zn, pn, A = (xn + λ) / 4, (yn + λ) / 4, (zn + λ) / 4, (pn + λ) / 4, (A + λ) / 4
        dxn, dyn, dzn, dpn, dA = (dxn + dλ) / 4, (dyn + dλ) / 4, (dzn + dλ) / 4, (dpn + dλ) / 4,
            (dA + dλ) / 4
        scale /= 4
    end
    X = (A0 - x) * scale / A
    Y = (A0 - y) * scale / A
    Z = (A0 - z) * scale / A
    dX = scale * (2 / 5 / A - (A0 - x) * dA / A^2)
    dY = scale * (2 / 5 / A - (A0 - y) * dA / A^2)
    dZ = scale * (2 / 5 / A - (A0 - z) * dA / A^2)
    P = -(X + Y + Z) / 2; dP = -(dX + dY + dZ) / 2
    XYZ = X * Y * Z; dXYZ = dX * Y * Z + X * dY * Z + X * Y * dZ
    E2 = X * Y + X * Z + Y * Z - 3 * P^2
    dE2 = dX * (Y + Z) + dY * (X + Z) + dZ * (X + Y) - 6 * P * dP
    E3 = XYZ + 2 * E2 * P + 4 * P^3
    dE3 = dXYZ + 2 * (dE2 * P + E2 * dP) + 12 * P^2 * dP
    E4 = (2 * XYZ + E2 * P + 3 * P^3) * P
    dE4 = (2 * dXYZ + dE2 * P + E2 * dP + 9 * P^2 * dP) * P + (2 * XYZ + E2 * P + 3 * P^3) * dP
    E5 = XYZ * P^2
    dE5 = dXYZ * P^2 + 2 * XYZ * P * dP
    series = 1 - 3 * E2 / 14 + E3 / 6 + 9 * E2^2 / 88 - 3 * E4 / 22 - 9 * E2 * E3 / 52 + 3 * E5 / 26
    dseries = -3 * dE2 / 14 + dE3 / 6 + 9 * E2 * dE2 / 44 - 3 * dE4 / 22 -
        9 * (dE2 * E3 + E2 * dE3) / 52 + 3 * dE5 / 26
    head = scale * series / (A * sqrt(A))
    dhead = scale * (dseries - 1.5 * series * dA / A) / (A * sqrt(A))
    return head + 6 * total, dhead + 6 * dtotal
end
_carlson_rj(x, y, z, p) = _carlson_rj_dp(x, y, z, p)[1]

# The Legendre forms at (sin φ, cos φ), 0 ≤ φ ≤ π/2, with the complement m1 = 1 − m: cos φ is
# supplied directly, not as 1 − sin²φ, which keeps their digits as φ → π/2
_ellip_f(s, c, m1) = s * _carlson_rf(c^2, c^2 + m1 * s^2, 1.0)
_ellip_e(s, c, m1) = s * _carlson_rf(c^2, c^2 + m1 * s^2, 1.0) -
    (1 - m1) * s^3 * _carlson_rd(c^2, c^2 + m1 * s^2, 1.0) / 3
# D(φ|m) = (F − E)/m without the cancellation of F − E
_ellip_d(s, c, m1) = s^3 * _carlson_rd(c^2, c^2 + m1 * s^2, 1.0) / 3
# (Π(n; φ|m) − F(φ|m))/n, the primitive of the pole terms, with n1 = 1 − n (principal value
# where the characteristic n > 1 passes its pole)
_ellip_pole(s, c, m1, n1) = s^3 * _carlson_rj(c^2, c^2 + m1 * s^2, 1.0, c^2 + n1 * s^2) / 3
_ellip_pi(s, c, m1, n, n1) = _ellip_f(s, c, m1) + n * _ellip_pole(s, c, m1, n1)
# the second-order integral ∫ dθ/((1 − n sin²θ)² √(1 − m sin²θ)) = Π + n ∂Π/∂n, the derivative
# from ∂R_J/∂p: no coefficient 1/(n − 1) or 1/(m − n) (the classical reduction loses every digit
# as a root or the point at infinity approaches the pole of the substitution, E → 1)
function _ellip_pi2(s, c, m1, n, n1)
    rj, drj = _carlson_rj_dp(c^2, c^2 + m1 * s^2, 1.0, c^2 + n1 * s^2)
    return _ellip_f(s, c, m1) + n * s^3 * (2 * rj - n * s^2 * drj) / 3
end

# the complete integrals
_ellip_k(m1) = _carlson_rf(0.0, m1, 1.0)
_ellip_e_complete(m1) = _carlson_rf(0.0, m1, 1.0) - (1 - m1) * _carlson_rd(0.0, m1, 1.0) / 3
_ellip_pi_complete(m1, n, n1) = _ellip_k(m1) + n * _carlson_rj(0.0, m1, 1.0, n1) / 3
function _ellip_pi2_complete(m1, n, n1)
    rj, drj = _carlson_rj_dp(0.0, m1, 1.0, n1)
    return _ellip_k(m1) + n * (2 * rj - n * drj) / 3
end

# the same forms for any real amplitude: φ = kπ + φ₀ with |φ₀| ≤ π/2 gives 2k × the complete
# integral plus the value at φ₀ (sin φ₀ signed, cos φ₀ ≥ 0). The complete integral enters only
# once a half period is crossed: for n = 1 it diverges while the incomplete value is finite.
function _legendre(term, complete, φ)
    k = round(φ / pi)
    s, c = sincos(φ - k * pi)
    return k == 0 ? term(s, c) : 2 * k * complete() + term(s, c)
end
_ellip_f(φ, m1) = _legendre((s, c) -> _ellip_f(s, c, m1), () -> _ellip_k(m1), φ)
_ellip_e(φ, m1) = _legendre((s, c) -> _ellip_e(s, c, m1), () -> _ellip_e_complete(m1), φ)
_ellip_d(φ, m1) = _legendre((s, c) -> _ellip_d(s, c, m1), () -> _carlson_rd(0.0, m1, 1.0) / 3, φ)
_ellip_pole(φ, m1, n1) = _legendre((s, c) -> _ellip_pole(s, c, m1, n1),
    () -> _carlson_rj(0.0, m1, 1.0, n1) / 3, φ)
_ellip_pi(φ, m1, n, n1) = _legendre((s, c) -> _ellip_pi(s, c, m1, n, n1),
    () -> _ellip_pi_complete(m1, n, n1), φ)
_ellip_pi2(φ, m1, n, n1) = _legendre((s, c) -> _ellip_pi2(s, c, m1, n, n1),
    () -> _ellip_pi2_complete(m1, n, n1), φ)

# the APEX interfaces pass the parameter itself; the complement is formed here
_complete_pi(n, n1, m) = _ellip_pi_complete(1 - m, n, n1)
_elliptic_D(φ, m) = _ellip_d(φ, 1 - m)

# ∫_r^∞ dr/√(lead (r − x1)(r − x2)(r − x3)(r − x4)) for four roots below r, real or a complex
# pair x3 = conj(x4), in Carlson's form 2 R_F(U12², U13², U14²)/√lead, U_ij = Y_i Y_j + Y_k Y_l,
# Y_i = √(r − x_i): the terms are positive (or complex conjugates), so the Mino time to infinity
# keeps its digits
function _mino_to_infinity(lead, roots, r)
    Y1, Y2, Y3, Y4 = sqrt(r - roots[1]), sqrt(r - roots[2]), sqrt(r - roots[3]), sqrt(r - roots[4])
    return 2 * real(_carlson_rf((Y1 * Y2 + Y3 * Y4)^2, (Y1 * Y3 + Y2 * Y4)^2,
        (Y1 * Y4 + Y2 * Y3)^2)) / sqrt(lead)
end

"""
    _landen(m, m1; complementary_modulus=sqrt(m1))

The descending Landen (arithmetic–geometric mean) sequence of the parameter `m`, started
from its complement `m1` = 1 − m as the caller formed it: a₀ = 1, b₀ = √m1. Returns
`(m, m1, a, b, K)` with K(m) = π/(2 AGM(1, √m1)). The caller can retain
`complementary_modulus` before its square underflows.
"""
function _landen(m, m1; complementary_modulus=sqrt(m1))
    a, b = 1.0, float(complementary_modulus)
    as, bs = Float64[], Float64[]
    for _ in 1:64
        push!(as, a); push!(bs, b)
        # the next mean is the limit to (a − b)²/(8a²) < eps
        abs(a - b) <= 1.0e-8 * a && return (m=m, m1=m1, a=as, b=bs, K=pi / (a + b))
        a, b = (a + b) / 2, sqrt(a * b)
    end
    error("The Landen sequence of m1 = $m1 does not converge.")
end

# (sn, cn, dn)(u | m) for 0 ≤ u ≤ K/2: Bulirsch's recursion for cot(am u) down the Landen
# sequence. Every step combines positive numbers, so all three keep relative accuracy.
function _sncndn(u, L)
    n = length(L.a)
    mean = (L.a[n] + L.b[n]) / 2
    s, c = sincos(u * mean)
    iszero(s) && return (0.0, 1.0, 1.0)
    t = c / s
    c = mean * t
    dn = 1.0
    for i in n:-1:1
        t *= c
        c *= dn
        dn = (L.b[i] + t) / (L.a[i] + t)
        t = c / L.a[i]
    end
    h = hypot(c, 1.0)
    return (1 / h, c / h, dn)
end

"""
(sn, cn, dn)(u | m) for the parameter record `L` of `_landen`. u is reduced to [−K, K] (sn,
cn have period 4K and change sign under u → u ± 2K); on (K/2, K] the quarter-period
reflection sn(K − w) = cd w, cn(K − w) = k' sd w, dn(K − w) = k' nd w keeps cn and dn
relatively accurate where they are small.
"""
function _ellipj_reduced(u, L)
    K = L.K
    y = u - 4K * round(u / (4K))                 # [−2K, 2K]
    flip = false
    if y > K
        y -= 2K; flip = true
    elseif y < -K
        y += 2K; flip = true
    end
    v = abs(y)
    sn, cn, dn = if v <= K / 2
        _sncndn(v, L)
    else
        sw, cw, dw = _sncndn(K - v, L)
        kp = L.b[1]
        (cw / dw, kp * sw / dw, kp / dw)
    end
    sn = copysign(sn, y)
    return flip ? (-sn, -cn, dn) : (sn, cn, dn)
end
