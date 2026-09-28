# Closed-form r(λ) and λ(r) shared by the plunges B1, B3, B4 and the Trapped cases N1–N6
# (t, φ, τ come from the radial engine). λ = 0 at the turning point `turn`.

function _four_simple_inner_model(a, energy, roots)
    x1, x2, x3, x4 = roots
    kappa = 1.0 - energy^2
    n = (x2 - x1) / (x3 - x1)
    m = (x4 - x3) * (x2 - x1) /
        ((x4 - x2) * (x3 - x1))
    xi = sqrt(kappa * (x4 - x2) * (x3 - x1)) / 2.0
    b = x2 - x3

    function radius(q)
        u = xi * q
        sn = Elliptic.Jacobi.sn(u, m)
        return x3 + b / (1.0 - n * sn^2)
    end
    function q_of_radius(r)
        s2 = (1.0 - b / (r - x3)) / n
        phi = asin(sqrt(clamp(s2, 0.0, 1.0)))
        return Elliptic.F(phi, m) / xi
    end
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:four_simple_inner,
        turn=x2,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        modulus=m,
    )
end

function _outer_double_model(a, energy, roots)
    x1, x2, x3 = roots
    kappa = 1.0 - energy^2
    b = x2 - x1
    c = x3 - x1
    n = b / c
    scale = 2.0 / (sqrt(kappa) * c)
    phi_of_r(r) = asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))

    i0(r) = scale * _pi_real(n, phi_of_r(r), 0.0)
    i0_turn = i0(x2)
    function radius(q)
        target = i0_turn - q
        k = sqrt(1.0 - n)
        angle = clamp(k * target / scale, 0.0, pi / 2.0)
        tangent = tan(angle) / k
        s2 = tangent^2 / (1.0 + tangent^2)
        return x1 + b * s2
    end
    q_of_radius(r) = i0_turn - i0(r)
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:outer_double,
        turn=x2,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, x3=x3),
        modulus=0.0,
    )
end

function _four_simple_single_exterior_model(a, energy, roots)
    x1, x2, x3, x4 = roots
    kappa = 1.0 - energy^2
    n = (x4 - x3) / (x4 - x2)
    m = (x4 - x3) * (x2 - x1) /
        ((x4 - x2) * (x3 - x1))
    xi = sqrt(kappa * (x4 - x2) * (x3 - x1)) / 2.0
    b = x3 - x2
    kcomplete = Elliptic.K(m)

    function radius(q)
        u = clamp(kcomplete - xi * q, 0.0, kcomplete)
        sn = Elliptic.Jacobi.sn(u, m)
        return x2 + b / (1.0 - n * sn^2)
    end
    function q_of_radius(r)
        s2 = (1.0 - b / (r - x2)) / n
        u = Elliptic.F(asin(sqrt(clamp(s2, 0.0, 1.0))), m)
        return (kcomplete - u) / xi
    end
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:four_simple_single_exterior,
        turn=x4,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        modulus=m,
    )
end

function _complex_pair_model(a, energy, roots, rho, eta)
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

    # r = [bb x2 (1 + cn) + aa x1 (1 − cn)] / [bb (1 + cn) + aa (1 − cn)]: every term keeps
    # its sign, so neither a far turning point (x2 ~ 1/(1 − E²)) nor a far complex pair
    # costs digits; 1 ± cn is formed from sn² on the side where it is small
    function radius(q)
        sn, cn, _ = Elliptic.ellipj(xi * q, m)
        onepc = cn >= 0 ? 1 + cn : sn^2 / (1 - cn)
        onemc = cn <= 0 ? 1 - cn : sn^2 / (1 + cn)
        return (bb * x2 * onepc + aa * x1 * onemc) / (bb * onepc + aa * onemc)
    end
    function q_of_radius(r)
        y = (bb * (x2 - r) - aa * (r - x1)) /
            (bb * (x2 - r) + aa * (r - x1))
        return Elliptic.F(pi / 2.0 + asin(clamp(y, -1.0, 1.0)), m) / xi
    end
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:two_real_complex_pair,
        turn=x2,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, rho=rho, eta=eta, A=aa, B=bb),
        modulus=m,
    )
end

function _parabolic_three_simple_model(a, roots)
    x1, x2, x3 = roots
    b = x2 - x1
    c = x3 - x1
    m = b / c
    scale = sqrt(2.0 / c)
    phi_of_r(r) = asin(sqrt(clamp((r - x1) / b, 0.0, 1.0)))

    i0(r) = scale * Elliptic.F(phi_of_r(r), m)
    i0_turn = i0(x2)
    function radius(q)
        target = clamp((i0_turn - q) / scale, 0.0, Elliptic.K(m))
        sn = Elliptic.Jacobi.sn(target, m)
        return x1 + b * sn^2
    end
    q_of_radius(r) = i0_turn - i0(r)
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:parabolic_three_simple,
        turn=x2,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, x3=x3),
        modulus=m,
    )
end

function _hyperbolic_four_simple_model(a, energy, roots)
    x1, x2, x3, x4 = roots
    n = (x3 - x2) / (x4 - x2)
    m = (x3 - x2) * (x4 - x1) /
        ((x4 - x2) * (x3 - x1))
    g = x3 - x4
    scale = 2.0 / sqrt(_e2m1(energy) *
        (x4 - x2) * (x3 - x1))
    function radius(q)
        u = clamp(q / scale, 0.0, Elliptic.K(m))
        sn = Elliptic.Jacobi.sn(u, m)
        return x4 + g / (1.0 - n * sn^2)
    end
    function q_of_radius(r)
        s2 = (x4 - x2) * (x3 - r) /
            ((x3 - x2) * (x4 - r))
        phi = asin(sqrt(clamp(s2, 0.0, 1.0)))
        return scale * Elliptic.F(phi, m)
    end
    lambda_horizon = q_of_radius(_rplus(a))
    return (
        kind=:hyperbolic_four_simple,
        turn=x3,
        lambda_horizon=lambda_horizon,
        r=radius,
        lambda_of_r=q_of_radius,
        roots=(x1=x1, x2=x2, x3=x3, x4=x4),
        modulus=m,
    )
end

function _simple_root_data(a, energy, lz, q)
    structure = kerr_geo_root_structure(a, energy, lz, q)
    return sort(Float64[root.radius for root in structure.real_roots]), _nonreal_roots(structure)
end

"""The radial model of B1, B3 or B4."""
function _plunge_radial_model(case_id, a, energy, lz, q)
    roots, nonreal = _simple_root_data(a, energy, lz, q)
    case_id === :B1 && return _four_simple_inner_model(
        a, energy, roots)
    case_id === :B3 && return _four_simple_single_exterior_model(
        a, energy, roots)
    if case_id === :B4
        length(nonreal) == 2 || error("B4 requires one conjugate radial-root pair.")
        return _complex_pair_model(a, energy, roots,
            real(first(nonreal)), abs(imag(first(nonreal))))
    end
    error("No shared radial model for $(case_id).")
end

"""The radial model of the Trapped case with disposition NFD0k (= case Nk)."""
function _trapped_radial_model(disposition, a, energy, lz, q)
    roots, nonreal = _simple_root_data(a, energy, lz, q)
    disposition === :NFD01 && return _four_simple_inner_model(a, energy, roots)
    disposition === :NFD02 && return _outer_double_model(a, energy, roots)
    disposition === :NFD03 && return _four_simple_single_exterior_model(
        a, energy, roots)
    if disposition === :NFD04
        length(nonreal) == 2 || error("NFD04 requires one conjugate root pair.")
        return _complex_pair_model(a, energy, roots,
            real(first(nonreal)), abs(imag(first(nonreal))))
    end
    disposition === :NFD05 && return _parabolic_three_simple_model(a, roots)
    disposition === :NFD06 && return _hyperbolic_four_simple_model(
        a, energy, roots)
    error("Unknown Trapped disposition $(disposition).")
end
