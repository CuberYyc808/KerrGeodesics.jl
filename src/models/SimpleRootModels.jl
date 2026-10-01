# Closed-form radial models of the plunges B1–B6 and the Trapped cases N1–N6 (t, φ, τ come from
# the radial engine). Every model gives `radius(δ)`, the radius at Mino time δ ≥ 0 after the
# turning point `turn` (the leg runs inward, `inward = true`), and its inverse `mino(r)`; B2, B5
# and B6 also carry `basis(r) = (I0, I1, I2)` and `pole(h, r)` for the analytic
# primitive checks and radial inverse maps.

# A1: a periodic libration, starting at the inner turning point. The same radial
# state is used by the ordinary and exact-extremal members.
function _libration_model(energy, r1, r2, r3, r4)
    m = (r1-r2)*(r3-r4)/((r1-r3)*(r2-r4))
    m1 = (r2-r3)*(r1-r4)/((r1-r3)*(r2-r4))
    L = _landen(m,m1)
    omega = sqrt(max(0.0,(1-energy)*(1+energy)*(r1-r3)*(r2-r4)))/2
    function state(lambda)
        sn,cn,dn = _ellipj_reduced(omega*lambda,L)
        D = (r1-r3)*cn^2+(r2-r3)*sn^2
        w = (r1-r2)*(r2-r3)
        return (r2+w*sn^2/D,2omega*w*(r1-r3)*sn*cn*dn/D^2)
    end
    function mino(r)
        s2 = (r1-r3)*(r-r2)/((r1-r2)*(r-r3))
        return _ellip_f(asin(sqrt(clamp(s2,0.0,1.0))),m1)/omega
    end
    return (kind=:four_real_libration, radius=lambda->state(lambda)[1],
        velocity=lambda->state(lambda)[2], state=state, mino=mino, inward=false,
        period=2L.K/omega, omega=omega, K=L.K)
end

# B1, N1: four simple roots, the inner allowed interval from x2 inward
function _four_simple_inner_model(energy, roots)
    x1, x2, x3, x4 = roots
    kappa = -_e2m1(energy)
    n = (x2 - x1) / (x3 - x1)
    m = (x4 - x3) * (x2 - x1) /
        ((x4 - x2) * (x3 - x1))
    m1 = (x3 - x2) * (x4 - x1) / ((x4 - x2) * (x3 - x1))
    L = _landen(m, m1)
    xi = sqrt(kappa * (x4 - x2) * (x3 - x1)) / 2.0
    b = x2 - x3
    function radius(δ)
        sn = _ellipj_reduced(xi * δ, L)[1]
        return x3 + b / (1.0 - n * sn^2)
    end
    function mino(r)
        s2 = (1.0 - b / (r - x3)) / n
        phi = asin(sqrt(clamp(s2, 0.0, 1.0)))
        return _ellip_f(phi, m1) / xi
    end
    return (kind=:four_simple_inner, turn=x2, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4), modulus=m)
end

# B2, N2: the double root x3 above the simple roots x1 < x2; motion between x1 and x2
function _b2_model(energy, roots)
    x1, x2, x3 = roots
    kappa = -_e2m1(energy)
    b = x2 - x1
    c = x3 - x1
    n = b / c
    n1 = (x3 - x2) / c
    scale = 2 / (sqrt(kappa) * c)
    phi_of_r(r) = asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))
    function basis(r)
        phi = phi_of_r(r)
        pin = _ellip_pi(phi, 1.0, n, n1)
        i0 = scale * pin
        i1 = scale * (x3 * pin - c * phi)
        i2 = scale * (x3^2 * pin - 2 * x3 * c * phi +
            c^2 * ((1 - n / 2) * phi + n * sin(2 * phi) / 4))
        return (I0=i0, I1=i1, I2=i2)
    end
    function pole(h, r)
        phi = phi_of_r(r)
        nh = b / (h - x1)
        return scale / (x1 - h) * (
            -n / (nh - n) * _ellip_pi(phi, 1.0, n, n1) +
            nh / (nh - n) * _ellip_pi(phi, 1.0, nh, (h - x2) / (h - x1)))
    end
    i0_turn = basis(x2).I0
    function radius(δ)
        k = sqrt(1 - n)
        angle = clamp(k * (i0_turn - δ) / scale, 0.0, pi / 2)
        tangent = tan(angle) / k
        s2 = tangent^2 / (1 + tangent^2)
        return x1 + b * s2
    end
    mino(r) = i0_turn - basis(r).I0
    return (kind=:b2_outer_double_inner_simple, turn=x2, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, x3=x3), modulus=0.0, basis=basis, pole=pole)
end

# B3, N3: four simple roots, the single exterior root x4 as the turning point
function _four_simple_single_exterior_model(energy, roots)
    x1, x2, x3, x4 = roots
    kappa = -_e2m1(energy)
    n = (x4 - x3) / (x4 - x2)
    m = (x4 - x3) * (x2 - x1) /
        ((x4 - x2) * (x3 - x1))
    m1 = (x3 - x2) * (x4 - x1) / ((x4 - x2) * (x3 - x1))
    L = _landen(m, m1)
    xi = sqrt(kappa * (x4 - x2) * (x3 - x1)) / 2.0
    function radius(δ)
        sn,cn,dn = _ellipj_reduced(clamp(xi*δ,0.0,L.K),L)
        # sn(K-u)=cn(u)/dn(u); both terms in the denominator are positive.
        return x2+(x4-x2)*dn^2/(cn^2+(x4-x1)/(x3-x1)*sn^2)
    end
    function mino(r)
        phi=atan(sqrt(max(0.0,(x4-r)*(x3-x1))),
            sqrt(max(0.0,(r-x3)*(x4-x1))))
        return _ellip_f(phi,m1)/xi
    end
    return (kind=:four_simple_single_exterior, turn=x4, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4), modulus=m)
end

# B4, N4: two real roots x1 < x2 and a complex pair ρ ± iη; the turning point x2
function _complex_pair_model(energy, roots, rho, eta)
    x1, x2 = roots
    aa = hypot(x2 - rho, eta)
    bb = hypot(x1 - rho, eta)
    # aa − bb from aa² − bb² = (x2 − x1)(x2 + x1 − 2ρ): no cancellation when the complex
    # pair is far away (ρ ~ 1/(1 − E²) as E → 1)
    dab = (x2 - x1) * (x2 + x1 - 2rho) / (aa + bb)
    xi = sqrt(-_e2m1(energy) * aa * bb)
    # m = ((x2 − x1)² − (aa − bb)²)/(4 aa bb), the two factors formed from
    # w(x) = |x − z| − (x − ρ) and w̄(x) = |x − z| + (x − ρ) without cancellation
    wminus(x, h) = x - rho > 0 ? eta^2 / (h + (x - rho)) : h - (x - rho)
    wplus(x, h) = x - rho < 0 ? eta^2 / (h - (x - rho)) : h + (x - rho)
    m = (wminus(x1, bb) - wminus(x2, aa)) * (wplus(x2, aa) - wplus(x1, bb)) / (4.0 * aa * bb)
    # 1 − m = ((aa + bb)² − (x2 − x1)²)/(4 aa bb), its first factor the sum of the two w's
    m1 = (wminus(x2, aa) + wplus(x1, bb)) * (aa + bb + (x2 - x1)) / (4.0 * aa * bb)
    L = _landen(m, m1)
    # r = [bb x2 (1 + cn) + aa x1 (1 − cn)] / [bb (1 + cn) + aa (1 − cn)]: every term keeps
    # its sign, so neither a far turning point (x2 ~ 1/(1 − E²)) nor a far complex pair
    # costs digits; 1 ± cn is formed from sn² on the side where it is small
    function radius(δ)
        sn, cn, _ = _ellipj_reduced(xi * δ, L)
        onepc = cn >= 0 ? 1 + cn : sn^2 / (1 - cn)
        onemc = cn <= 0 ? 1 - cn : sn^2 / (1 + cn)
        return (bb * x2 * onepc + aa * x1 * onemc) / (bb * onepc + aa * onemc)
    end
    function mino(r)
        y = (bb * (x2 - r) - aa * (r - x1)) /
            (bb * (x2 - r) + aa * (r - x1))
        return _ellip_f(pi / 2.0 + asin(clamp(y, -1.0, 1.0)), m1) / xi
    end
    return (kind=:two_real_complex_pair, turn=x2, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, rho=rho, eta=eta, A=aa, B=bb), modulus=m)
end

# B5, N5: E = 1, three simple roots; motion between x1 and x2
function _b5_model(roots)
    x1, x2, x3 = roots
    b = x2 - x1
    c = x3 - x1
    m = b / c
    m1 = (x3 - x2) / c
    L = _landen(m, m1)
    scale = sqrt(2 / c)
    phi_of_r(r) = asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))
    j2(phi) = _ellip_d(phi, m1)                     # (F − E)/m
    function j4(phi)
        boundary = sin(phi) * cos(phi) *
            sqrt(max(1 - m * sin(phi)^2, 0.0))
        return (boundary - _ellip_f(phi, m1) +
            (2 + 2 * m) * j2(phi)) / (3 * m)
    end
    function basis(r)
        phi = phi_of_r(r)
        f = _ellip_f(phi, m1)
        return (
            I0=scale * f,
            I1=scale * (x1 * f + b * j2(phi)),
            I2=scale * (x1^2 * f + 2 * x1 * b * j2(phi) + b^2 * j4(phi)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        return scale * _ellip_pi(phi, m1, b / (h - x1), (h - x2) / (h - x1)) / (x1 - h)
    end
    i0_turn = basis(x2).I0
    function radius(δ)
        sn = _ellipj_reduced(clamp((i0_turn - δ) / scale, 0.0, L.K), L)[1]
        return x1 + b * sn^2
    end
    mino(r) = i0_turn - basis(r).I0
    return (kind=:b5_parabolic_three_simple, turn=x2, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, x3=x3), modulus=m, basis=basis, pole=pole)
end

# B6, N6: E > 1, four simple roots; motion between x2 and x3
function _b6_model(energy, roots)
    x1, x2, x3, x4 = roots
    lead = _e2m1(energy)
    n = (x3 - x2) / (x4 - x2)
    n1 = (x4 - x3) / (x4 - x2)
    m = (x3 - x2) * (x4 - x1) / ((x4 - x2) * (x3 - x1))
    m1 = (x2 - x1) * (x4 - x3) / ((x4 - x2) * (x3 - x1))
    L = _landen(m, m1)
    g = x3 - x4
    scale = 2 / sqrt(lead * (x4 - x2) * (x3 - x1))
    function phi_of_r(r)
        s2 = (x4 - x2) * (x3 - r) / ((x3 - x2) * (x4 - r))
        return asin(sqrt(clamp(s2, 0.0, 1.0)))
    end
    function basis(r)
        phi = phi_of_r(r)
        f = _ellip_f(phi, m1)
        pin = _ellip_pi(phi, m1, n, n1)
        return (
            I0=-scale * f,
            I1=-scale * (x4 * f + g * pin),
            I2=-scale * (x4^2 * f + 2 * x4 * g * pin +
                g^2 * _ellip_pi2(phi, m1, n, n1)),
        )
    end
    function pole(h, r)
        phi = phi_of_r(r)
        nh = n * (x4 - h) / (x3 - h)
        return -scale * (_ellip_f(phi, m1) / (x4 - h) + (x4 - x3) *
            _ellip_pi(phi, m1, nh, (x2 - h) * (x4 - x3) / ((x3 - h) * (x4 - x2))) /
            ((x4 - h) * (x3 - h)))
    end
    function radius(δ)
        s2 = _ellipj_reduced(clamp(δ / scale, 0.0, L.K), L)[1]^2
        numerator = (x4 - x2) * x3 - s2 * (x3 - x2) * x4
        denominator = (x4 - x2) - s2 * (x3 - x2)
        return numerator / denominator
    end
    mino(r) = scale * _ellip_f(phi_of_r(r), m1)
    return (kind=:b6_hyperbolic_four_simple, turn=x3, radius=radius, mino=mino, inward=true,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4), modulus=m, basis=basis, pole=pole)
end

_simple_root_data(structure) =
    (sort(Float64[root.radius for root in structure.real_roots]), _nonreal_roots(structure))

"""The radial model of B1–B6 from the classified root structure."""
function _plunge_radial_model(case_id, energy, structure)
    roots, nonreal = _simple_root_data(structure)
    case_id === :B1 && return _four_simple_inner_model(energy, roots)
    case_id === :B2 && return _b2_model(energy, (roots[1], roots[2], roots[3]))
    case_id === :B3 && return _four_simple_single_exterior_model(energy, roots)
    if case_id === :B4
        length(nonreal) == 2 || error("B4 requires one conjugate radial-root pair.")
        return _complex_pair_model(energy, roots,
            real(first(nonreal)), abs(imag(first(nonreal))))
    end
    case_id === :B5 && return _b5_model((roots[1], roots[2], roots[3]))
    case_id === :B6 && return _b6_model(energy, (roots[1], roots[2], roots[3], roots[4]))
    error("No closed-form radial model for $(case_id).")
end

"""The radial model of the Trapped case with disposition NFD0k (= case Nk)."""
function _trapped_radial_model(disposition, energy, structure)
    roots, nonreal = _simple_root_data(structure)
    disposition === :NFD01 && return _four_simple_inner_model(energy, roots)
    disposition === :NFD02 && return _b2_model(energy, (roots[1], roots[2], roots[3]))
    disposition === :NFD03 && return _four_simple_single_exterior_model(energy, roots)
    if disposition === :NFD04
        length(nonreal) == 2 || error("NFD04 requires one conjugate root pair.")
        return _complex_pair_model(energy, roots,
            real(first(nonreal)), abs(imag(first(nonreal))))
    end
    disposition === :NFD05 && return _b5_model((roots[1], roots[2], roots[3]))
    disposition === :NFD06 && return _b6_model(energy, (roots[1], roots[2], roots[3], roots[4]))
    error("Unknown Trapped disposition $(disposition).")
end
