# Jacobi elliptic functions: the descending Landen (AGM) record of a parameter, Bulirsch's
# recursion for sn, cn, dn with the quarter-period reflection, the amplitude, and any real
# parameter through Abramowitz & Stegun 16.10 (m < 0) and 16.11 (m > 1). Every record is
# immutable once built, so one record serves concurrent evaluations.

"""
    _LandenRecord{T}

The descending Landen (arithmetic–geometric mean) sequence of the parameter `m`, started from
its complement `m1` = 1 − m as the caller formed it: a₀ = 1, b₀ = √m1, with
K(m) = π/(2 AGM(1, √m1)). For m1 = 0 the mean vanishes and K = ∞ (sn = tanh, cn = dn = sech);
a NaN parameter gives K = NaN and NaN functions.
"""
struct _LandenRecord{T}
    m::T
    m1::T
    a::Vector{T}
    b::Vector{T}
    K::T
end

"""
    _landen(m, m1; complementary_modulus=sqrt(m1))

The `_LandenRecord` of the parameter `m` from its complement `m1`. The caller can retain
`complementary_modulus` before its square underflows.
"""
function _landen(m, m1; complementary_modulus=sqrt(m1))
    T = _float_type(m, m1, complementary_modulus)
    a, b = one(T), T(complementary_modulus)
    as, bs = T[], T[]
    iszero(b) && return _LandenRecord{T}(m, m1, [a], [b], T(Inf))
    isnan(b) && return _LandenRecord{T}(m, m1, [a], [b], b)      # NaN propagates
    tol = _tol(T, 1.0e-8)
    for _ in 1:64
        push!(as, a); push!(bs, b)
        # the next mean is the limit to (a − b)²/(8a²) < eps
        abs(a - b) <= tol * a && return _LandenRecord{T}(m, m1, as, bs, pi / (a + b))
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
    iszero(s) && return (zero(s), one(s), one(s))
    t = c / s
    c = mean * t
    dn = one(t)
    for i in n:-1:1
        t *= c
        c *= dn
        dn = (L.b[i] + t) / (L.a[i] + t)
        t = c / L.a[i]
    end
    h = hypot(c, one(c))
    return (1 / h, c / h, dn)
end

"""
(sn, cn, dn)(u | m) for the parameter record `L` of `_landen`. u is reduced to [−K, K] (sn,
cn have period 4K and change sign under u → u ± 2K); on (K/2, K] the quarter-period
reflection sn(K − w) = cd w, cn(K − w) = k' sd w, dn(K − w) = k' nd w keeps cn and dn
relatively accurate where they are small. For m = 1 (K = ∞) they are tanh u, sech u, sech u.
"""
function _ellipj_reduced(u, L)
    K = L.K
    isinf(K) && return (tanh(u), sech(u), sech(u))
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

# the amplitude am(u | m), 0 ≤ m ≤ 1, continuous and increasing: am(u + 2jK) = am(u) + jπ, and
# on |y| ≤ K it is atan(sn y, cn y) ∈ [−π/2, π/2]
function _jacobi_am(u, L::_LandenRecord)
    iszero(u) && return zero(float(u))              # for every m, also a NaN parameter
    j = round(u / (2 * L.K))
    y = iszero(j) ? u : u - 2j * L.K
    sn, cn, _ = _ellipj_reduced(y, L)
    return j * pi + atan(sn, cn)
end

"""
    _JacobiParameter{T}

sn, cn, dn for any real parameter m from the Landen record of a parameter μ ∈ [0, 1]:
μ = m for 0 ≤ m ≤ 1; for m < 0 (A&S 16.10) μ = −m/(1 − m), and sn(u|m) = √μ1 sd(v|μ),
cn(u|m) = cd(v|μ), dn(u|m) = nd(v|μ) with v = u√(1 − m), μ1 = 1/(1 − m); for m > 1 (A&S
16.11) μ = 1/m, and sn(u|m) = √μ sn(v|μ), cn(u|m) = dn(v|μ), dn(u|m) = cn(v|μ) with v = u√m.
"""
struct _JacobiParameter{T}
    m::T
    scale::T            # v = scale × u
    factor::T           # √μ1 (m < 0), √μ (m > 1), 1 otherwise
    landen::_LandenRecord{T}
end

function _jacobi_parameter(m)
    T = float(typeof(m))
    if m < 0
        μ1 = 1 / (1 - m)
        return _JacobiParameter{T}(m, sqrt(1 - m), sqrt(μ1), _landen(-m * μ1, μ1))
    elseif m > 1
        μ = 1 / m
        return _JacobiParameter{T}(m, sqrt(m), sqrt(μ), _landen(μ, (m - 1) / m))
    end
    return _JacobiParameter{T}(m, one(T), one(T), _landen(m, 1 - m))
end

function _jacobi_sncndn(u, J::_JacobiParameter)
    iszero(u) && return (zero(float(u)), one(float(u)), one(float(u)))   # every m, also NaN
    sn, cn, dn = _ellipj_reduced(J.scale * u, J.landen)
    J.m < 0 && return (J.factor * sn / dn, cn / dn, 1 / dn)
    J.m > 1 && return (J.factor * sn, dn, cn)
    return (sn, cn, dn)
end
