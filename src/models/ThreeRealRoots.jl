# Legendre form of the parabolic (E = 1) radial motion outside three real roots x1 < x2 < x3,
# R = 2(r − x1)(r − x2)(r − x3) (the scattering leg of D1 and the capture leg of C2): with
# A = x3 − x1, m = (x2 − x1)/A and sin²φ = A/(r − x1), dλ = −√(2/A) dφ/√(1 − m sin²φ), so the
# Mino time from the infinity endpoint down to r is √(2/A) F(φ|m) and the leg reaches the
# turning point at λ∞ = √(2/A) K(m).

function _three_real_leg(roots)
    x1, x2, x3 = roots.x1, roots.x2, roots.x3
    A = x3 - x1
    m = (x2 - x1) / A
    m1 = (x3 - x2) / A
    L = _landen(m, m1)
    scale = sqrt(2 / A)
    return (x1=x1, x2=x2, x3=x3, A=A, m=m, m1=m1, L=L, scale=scale,
        lambda_infinity=scale * L.K)
end

# φ from tan²φ = A/(r − x3): π/2 at the turning point without asin next to 1
function _three_real_amplitude(leg, r)
    r >= leg.x3 - _radius_tol(_float_type(r, leg.x3)) || error("The radius lies below the turning point x3.")
    return atan(sqrt(leg.A / max(r - leg.x3, 0.0)))
end

# Mino time from the infinity endpoint down to r
_three_real_mino_from_infinity(leg, r) =
    leg.scale * _ellip_f(_three_real_amplitude(leg, r), leg.m1)

# r at Mino time δ after the infinity endpoint: r − x3 = A cn²u/sn²u with u = δ/scale is exact
# at the turning point (u = K), and r ≈ x1 + 2/δ² as δ → 0
function _three_real_radius_from_infinity(leg, δ)
    sn, cn, _ = _ellipj_reduced(δ / leg.scale, leg.L)
    return leg.x3 + leg.A * (cn / sn)^2
end

# the pole primitive (Π(n; φ|m) − F(φ|m))/(A n) of 1/(r − h), n = (h − x1)/A
_three_real_pole(leg, h, phi) = leg.scale * _ellip_pole(phi, leg.m1, (leg.x3 - h) / leg.A) / leg.A
