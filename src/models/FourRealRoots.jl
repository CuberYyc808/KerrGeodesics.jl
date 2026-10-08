# Legendre form of the radial motion outside four real roots r_A < r_B < r_C < r_D (the
# scattering leg of D2 and the capture leg of C4): modulus, amplitude, Mino-time prefactor, the
# Mino time to infinity and r from the Mino time measured from the infinity endpoint.

# the leg's constants: m and its complement from root differences, the Landen sequence, the
# prefactor, and sn²u∞ = s∞ with (sn, cn, dn)(u∞) in closed form
function _four_real_leg(energy, roots)
    rA, rB, rC, rD = roots.rA, roots.rB, roots.rC, roots.rD
    m = ((rD - rA) * (rC - rB)) / ((rD - rB) * (rC - rA))
    m1 = ((rB - rA) * (rD - rC)) / ((rD - rB) * (rC - rA))
    L = _landen(m, m1)
    prefactor = 2 / sqrt(_e2m1(energy) * (rD - rB) * (rC - rA))
    s∞ = (rC - rA) / (rD - rA)
    c∞2 = (rD - rC) / (rD - rA)
    jac∞ = (sqrt(s∞), sqrt(c∞2), sqrt(c∞2 + m1 * s∞))
    phi∞ = atan(sqrt((rC - rA) / (rD - rC)))
    u∞ = _ellip_f(phi∞, m1)
    return (roots=roots, m=m, m1=m1, L=L, prefactor=prefactor, s∞=s∞, jac∞=jac∞,
        u∞=u∞, lambda_infinity=prefactor * u∞)
end

# am(u) with sin² am = s2(r): as atan, since s2 → 1 − (r_D − r_C)/(r_D − r_A) at infinity (and
# asin near 1 loses half the digits when r_A is far away, E → 1⁺)
function _four_real_amplitude(leg, r)
    roots = leg.roots
    r >= roots.rD - _radius_tol(_float_type(r, roots.rD)) || error("The radius lies below the outer turning point r_D.")
    return atan(sqrt(max((roots.rC - roots.rA) * (r - roots.rD), 0.0) /
        ((roots.rD - roots.rC) * (r - roots.rA))))
end

# Mino time from the turning point r_D out to r
_four_real_mino_from_turn(leg, r) = leg.prefactor * _ellip_f(_four_real_amplitude(leg, r), leg.m1)

"""
r on the branch r ≥ r_D at u = F(φ|m) from the turning point and du = u∞ − u from infinity.
sn²u = s2(r) (a Möbius map) gives r = r_D + (r_D − r_C) sn²u/(s∞ − sn²u), s∞ = s2(∞) = sn²u∞;
the difference s∞ − sn²u = sn(u∞ + u) sn(u∞ − u) (1 − m sn²u s∞) is formed from du itself, so r
keeps its digits up to the conditioning of λ ↦ r as r → ∞. sn(u∞ + u) comes from the addition
formula (every term positive on [0, K]).
"""
function _four_real_radius(leg, u, du)
    roots = leg.roots
    sn, cn, dn = _ellipj_reduced(u, leg.L)
    s, c, d = leg.jac∞
    sn_sum = (s * cn * dn + sn * c * d) / (1 - leg.m * leg.s∞ * sn^2)
    gap = sn_sum * _ellipj_reduced(du, leg.L)[1] * (1 - leg.m * sn^2 * leg.s∞)
    return roots.rD + (roots.rD - roots.rC) * sn^2 / gap
end

# r at Mino time δ after the infinity endpoint
_four_real_radius_from_infinity(leg, δ) =
    _four_real_radius(leg, leg.u∞ - δ / leg.prefactor, δ / leg.prefactor)

# the characteristic of the pole of 1/(r − h) and its complement
_four_real_characteristic(leg, h) = (roots = leg.roots;
    ((roots.rD - roots.rA) * (roots.rC - h) / ((roots.rC - roots.rA) * (roots.rD - h)),
     (h - roots.rA) * (roots.rD - roots.rC) / ((roots.rC - roots.rA) * (roots.rD - h))))

# Every radial model offers `radius(δ)` (δ ≥ 0 from its anchor) and its inverse `mino(r)`; the
# exact-extremal branches step from a radius `left` by an increment `delta` of I0, which grows
# with r, so `mino` grows with −delta on an inward leg
function _basis_wrapper(name, basis, pole, radius, mino; inward=true)
    inverse_from(left, delta) = radius(mino(left) + (inward ? -delta : delta))
    return (name=name, basis=basis, pole=pole, inverse_from=inverse_from, mino=mino,
        horizon_root=false)
end
function _outer_four_real_model(energy, ascending_roots)
    x1, x2, x3, x4 = ascending_roots
    r1, r2, r3, r4 = x4, x3, x2, x1
    kappa = -_e2m1(energy)
    m = (r1-r2)*(r3-r4)/((r1-r3)*(r2-r4))
    m1 = (r2-r3)*(r1-r4)/((r1-r3)*(r2-r4))
    L = _landen(m,m1)
    n = (r1-r2)/(r1-r3)
    n1 = (r2-r3)/(r1-r3)
    d = r2-r3
    scale = 2/sqrt(kappa*(r1-r3)*(r2-r4))
    phi_of_r(r) = asin(sqrt(clamp(
        (r1-r3)*(r-r2)/((r1-r2)*(r-r3)), 0.0, 1.0)))
    function basis(r)
        phi=phi_of_r(r)
        f=_ellip_f(phi,m1)
        pin=_ellip_pi(phi,m1,n,n1)
        j2=_ellip_pi2(phi,m1,n,n1)
        return (I0=scale*f,
            I1=scale*(r3*f+d*pin),
            I2=scale*(r3^2*f+2r3*d*pin+d^2*j2))
    end
    function pole(h,r)
        phi=phi_of_r(r)
        nh=n*(r3-h)/(r2-h)
        f=_ellip_f(phi,m1)
        return scale*(f/(r3-h)-
            d*_ellip_pi(phi,m1,nh,(r1-h)*(r2-r3)/((r1-r3)*(r2-h)))/((r2-h)*(r3-h)))
    end
    # δ from the lower turning point r2 = x3 outward (I0 = 0 there)
    function radius(δ)
        sn=_ellipj_reduced(clamp(δ/scale,0.0,L.K),L)[1]
        s2=sn^2
        return (r2-r3*n*s2)/(1-n*s2)
    end
    mino(r)=basis(r).I0
    return _basis_wrapper(:extremal_four_real_outer,basis,pole,radius,mino;inward=false)
end
