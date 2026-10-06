# The Legendre forms F, E, D = (F − E)/m, Π, the pole primitive (Π − F)/n and the second-order
# integral Π₂ built on Carlson's integrals, their complete values, the reduction of any real
# amplitude, and the Mino time to infinity in Carlson form. Every form takes the complement
# m1 = 1 − m (and n1 = 1 − n for Π) as the caller forms it from its roots, so a parameter
# within rounding of 1 keeps its digits.

# The Legendre forms at (sin φ, cos φ), 0 ≤ φ ≤ π/2, with the complement m1 = 1 − m: cos φ is
# supplied directly, not as 1 − sin²φ, which keeps their digits as φ → π/2
_ellip_f(s, c, m1) = s * _carlson_rf(c^2, c^2 + m1 * s^2, one(m1))
# E(φ|m) as the sum of positive terms (DLMF 19.25.9, 0 ≤ m ≤ 1): the form F − m D cancels as
# m → 1, where E → sin φ while F and m D diverge
function _ellip_e(s, c, m1)
    Δ2 = c^2 + m1 * s^2
    return m1 * s * _carlson_rf(c^2, Δ2, one(m1)) +
        (1 - m1) * m1 * s^3 * _carlson_rd(c^2, one(m1), Δ2) / 3 + (1 - m1) * s * c / sqrt(Δ2)
end
# D(φ|m) = (F − E)/m without the cancellation of F − E
_ellip_d(s, c, m1) = s^3 * _carlson_rd(c^2, c^2 + m1 * s^2, one(m1)) / 3
# (Π(n; φ|m) − F(φ|m))/n, the primitive of the pole terms, with n1 = 1 − n (principal value
# where the characteristic n > 1 passes its pole)
_ellip_pole(s, c, m1, n1) = s^3 * _carlson_rj(c^2, c^2 + m1 * s^2, one(m1), c^2 + n1 * s^2) / 3
_ellip_pi(s, c, m1, n, n1) = _ellip_f(s, c, m1) + n * _ellip_pole(s, c, m1, n1)
# the second-order integral ∫ dθ/((1 − n sin²θ)² √(1 − m sin²θ)) = Π + n ∂Π/∂n, the derivative
# from ∂R_J/∂p: no coefficient 1/(n − 1) or 1/(m − n) (the classical reduction loses every digit
# as a root or the point at infinity approaches the pole of the substitution, E → 1)
function _ellip_pi2(s, c, m1, n, n1)
    rj, drj = _carlson_rj_dp(c^2, c^2 + m1 * s^2, one(m1), c^2 + n1 * s^2)
    return _ellip_f(s, c, m1) + n * s^3 * (2 * rj - n * s^2 * drj) / 3
end

# the complete integrals
_ellip_k(m1) = _carlson_rf(zero(m1), m1, one(m1))
# E(m) = 2 R_G(0, m1, 1) as a sum of positive terms: m1 R_F + m m1 R_D(0, 1, m1)/3 for m ≥ 0
# (DLMF 19.21.10 with the middle argument last), R_F − m R_D(0, m1, 1)/3 for m < 0
_ellip_e_complete(m1) = m1 <= 1 ?
    m1 * _carlson_rf(zero(m1), m1, one(m1)) + (1 - m1) * m1 * _carlson_rd(zero(m1), one(m1), m1) / 3 :
    _carlson_rf(zero(m1), m1, one(m1)) - (1 - m1) * _carlson_rd(zero(m1), m1, one(m1)) / 3
_ellip_d_complete(m1) = _carlson_rd(zero(m1), m1, one(m1)) / 3
_ellip_pi_complete(m1, n, n1) = _ellip_k(m1) + n * _carlson_rj(zero(m1), m1, one(m1), n1) / 3
function _ellip_pi2_complete(m1, n, n1)
    rj, drj = _carlson_rj_dp(zero(m1), m1, one(m1), n1)
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
_ellip_d(φ, m1) = _legendre((s, c) -> _ellip_d(s, c, m1), () -> _ellip_d_complete(m1), φ)
_ellip_pole(φ, m1, n1) = _legendre((s, c) -> _ellip_pole(s, c, m1, n1),
    () -> _carlson_rj(zero(m1), m1, one(m1), n1) / 3, φ)
_ellip_pi(φ, m1, n, n1) = _legendre((s, c) -> _ellip_pi(s, c, m1, n, n1),
    () -> _ellip_pi_complete(m1, n, n1), φ)
_ellip_pi2(φ, m1, n, n1) = _legendre((s, c) -> _ellip_pi2(s, c, m1, n, n1),
    () -> _ellip_pi2_complete(m1, n, n1), φ)

# ∫_r^∞ dr/√(lead (r − x1)(r − x2)(r − x3)(r − x4)) for four roots below r, real or a complex
# pair x3 = conj(x4), in Carlson's form 2 R_F(U12², U13², U14²)/√lead, U_ij = Y_i Y_j + Y_k Y_l,
# Y_i = √(r − x_i): the terms are positive (or complex conjugates), so the Mino time to infinity
# keeps its digits
function _mino_to_infinity(lead, roots, r)
    Y1, Y2, Y3, Y4 = sqrt(r - roots[1]), sqrt(r - roots[2]), sqrt(r - roots[3]), sqrt(r - roots[4])
    return 2 * real(_carlson_rf((Y1 * Y2 + Y3 * Y4)^2, (Y1 * Y3 + Y2 * Y4)^2,
        (Y1 * Y4 + Y2 * Y3)^2)) / sqrt(lead)
end
