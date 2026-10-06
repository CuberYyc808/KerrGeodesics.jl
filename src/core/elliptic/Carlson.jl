# Carlson's symmetric integrals R_F, R_C (with ∂R_C/∂y), R_D and R_J (with ∂R_J/∂p) by
# duplication, in the floating-point type of their arguments: the duplication stops when the
# truncated series is exact to eps(T), so the result has the precision of the inputs. R_J with
# a negative fourth argument is the Cauchy principal value.

# Carlson's R_F(x, y, z) by duplication (x, y, z ≥ 0, at most one of them zero, or x ≥ 0 and
# y, z complex conjugates); with two arguments zero the integral diverges
function _carlson_rf(x, y, z)
    T = _real_type(x, y, z)
    iszero(x) + iszero(y) + iszero(z) >= 2 && return T(Inf)
    A0 = (x + y + z) / 3
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (3eps(T))^(one(T) / 6)
    xn, yn, zn, A, scale = x, y, z, A0, one(T)
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
# gives both (the derivative's closed form (R_C − √x/y)/(2(x − y)) cancels there); |ε| < 1/4,
# so the terms fall below eps(T) within a number of terms proportional to the precision.
function _carlson_rc_dy(x, y)
    T = _real_type(x, y)
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
        value = zero(T); derivative = zero(T); term = one(T); dterm = -one(T)
        small = _tol(T, 1e-18)
        for k in 0:_nterms(T, 33) - 1
            value += term / (2k + 1)
            k > 0 && (derivative += k * dterm / (2k + 1))
            k > 0 && abs(term) < small && break
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
    T = _real_type(x, y, z)
    A0 = (x + y + 3z) / 5
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z)) / (eps(T) / 4)^(one(T) / 6)
    xn, yn, zn, A, scale, total = x, y, z, A0, one(T), zero(T)
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
    T = _real_type(x, y, z, p)
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
    Q = max(abs(A0 - x), abs(A0 - y), abs(A0 - z), abs(A0 - p)) / (eps(T) / 4)^(one(T) / 6)
    xn, yn, zn, pn, A = x, y, z, p, A0
    dxn = dyn = dzn = zero(T); dpn = one(T); dA = T(2) / 5
    scale, total, dtotal = one(T), zero(T), zero(T)
    while scale * Q > abs(A)
        sx, sy, sz, sp = sqrt(xn), sqrt(yn), sqrt(zn), sqrt(pn)
        # an argument that is zero (a complete integral) has no p dependence at that step
        dsx, dsy, dsz = (iszero(dxn) ? zero(T) : dxn / (2sx)), (iszero(dyn) ? zero(T) : dyn / (2sy)),
            (iszero(dzn) ? zero(T) : dzn / (2sz))
        dsp = dpn / (2sp)
        λ = sx * sy + sx * sz + sy * sz
        dλ = dsx * (sy + sz) + dsy * (sx + sz) + dsz * (sx + sy)
        d = (sp + sx) * (sp + sy) * (sp + sz)
        dd = (dsp + dsx) * (sp + sy) * (sp + sz) + (sp + sx) * (dsp + dsy) * (sp + sz) +
            (sp + sx) * (sp + sy) * (dsp + dsz)
        arg = 1 + scale^3 * δ / d^2
        darg = scale^3 * (dδ / d^2 - 2 * δ * dd / d^3)
        rc, drc = _carlson_rc_dy(one(T), arg)
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
    dX = scale * (T(2) / 5 / A - (A0 - x) * dA / A^2)
    dY = scale * (T(2) / 5 / A - (A0 - y) * dA / A^2)
    dZ = scale * (T(2) / 5 / A - (A0 - z) * dA / A^2)
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
