# Radial models of the Critical members: r(λ) and the Mino-time primitives I0, I1, I2 of the
# motion that tends to a repeated root r_c (K2: triple root; K4/K5: E < 1 double root, the
# homoclinic substitution; K7/K8: E = 1 double root; K10/K11: E > 1 double root).

"""
Radial model of the Critical case `case_id` from the real roots of R (increasing; for K4/K5
the factorization (x1, r_c, r_a) of R = (1 − E²)(r − x1)(r − r_c)²(r_a − r)). Every model has
`basis(r) = (I0, I1, I2)`, `pole(h, r)`, `inverse_i0` and its radial range `lower`/`upper`.
"""
function _critical_radial_model(case_id, energy, radii)
    case_id === :K2 && return _k2_model(energy, (radii[1], radii[2]))
    case_id === :K8 && return _k8_model((radii[1], radii[2]))
    case_id === :K11 && return _k11_model(energy, (radii[1], radii[2], radii[3]))
    case_id === :K7 && return _k7_model((radii[1], radii[2]))
    case_id === :K10 && return _k10_model(energy, (radii[1], radii[2], radii[3]))
    h = _homoclinic_model(energy, radii...)
    # K4: the outer branch (r_c, r_a]; K5: the inner branch (x1, r_c), I0 increasing inwards
    case_id === :K4 && return (kind=:homoclinic_outer_branch, roots=h.roots, lower=h.lower,
        upper=h.upper, basis=h.basis, pole=h.pole, inverse_i0=h.inverse_i0)
    case_id === :K5 && return (kind=:homoclinic_inner_branch, roots=h.roots, upper=h.lower,
        basis=h.inner_basis, pole=h.inner_pole, inverse_i0=h.inner_inverse_i0)
    error("No Critical radial model is registered for $(case_id).")
end

function _k2_model(energy, roots)
    x1, x3 = roots
    kappa = -_e2m1(energy)
    d = x3 - x1
    scale = 2 / (d * sqrt(kappa))
    y_of_r(r) = sqrt(max((r - x1) / (x3 - r), 0.0))
    function basis(r)
        y = y_of_r(r)
        return (
            I0=scale * y,
            I1=scale * (x3 * y - d * atan(y)),
            I2=scale * (x3^2 * y - 2 * x3 * d * atan(y) +
                d^2 * (atan(y) + y / (1 + y^2)) / 2),
        )
    end
    function pole(h, r)
        y = y_of_r(r)
        return scale * (y / (x3 - h) + d / (x3 - h) *
            _quadratic_denominator_primitive(x1 - h, x3 - h, y))
    end
    function inverse_i0(target)
        y = max(target / scale, 0.0)
        return (x1 + x3 * y^2) / (1 + y^2)
    end
    return (kind=:k2_outer_triple, roots=(x1=x1, x3=x3), upper=x3,
        basis=basis, pole=pole, inverse_i0=inverse_i0)
end

function _k8_model(roots)
    x1, x2 = roots
    d = x2 - x1
    scale = sqrt(2 / d)
    y_of_r(r) = sqrt(clamp((r - x1) / d, 0.0, 1.0))
    function basis(r)
        y = y_of_r(r)
        at = _real_atanh(y)
        return (
            I0=scale * at,
            I1=scale * (x2 * at - d * y),
            I2=scale * (x2^2 * at - 2 * x2 * d * y +
                d^2 * (y - y^3 / 3)),
        )
    end
    function pole(h, r)
        y = y_of_r(r)
        at = _real_atanh(y)
        return scale * (at / (x2 - h) + d / (x2 - h) *
            _quadratic_denominator_primitive(x1 - h, d, y))
    end
    inverse_i0(target) = begin
        y = min(tanh(max(target / scale, 0.0)), prevfloat(1.0))
        x1 + d * y^2
    end
    return (kind=:k8_parabolic_outer_double, roots=(x1=x1, x2=x2), upper=x2,
        basis=basis, pole=pole, inverse_i0=inverse_i0)
end

function _k11_model(energy, roots)
    x1, x2, x3 = roots
    lead = _e2m1(energy)
    d = x2 - x1
    b = x3 - x2
    c = x3 - x1
    n = c / b
    scale = 2 / (sqrt(lead) * b)
    y_of_r(r) = sqrt(max((r - x2) / (r - x1), 0.0))
    j1sq(y) = y / (2 * (1 - y^2)) + _real_atanh(y) / 2
    function prod(y)
        return _j_inv(1.0, y) / (1 - n) - n * _j_inv(n, y) / (1 - n)
    end
    function second(y)
        return -n * _j_inv(1.0, y) / (n - 1)^2 -
            j1sq(y) / (n - 1) + n^2 * _j_inv(n, y) / (n - 1)^2
    end
    function basis(r)
        y = y_of_r(r)
        jn = _j_inv(n, y)
        return (
            I0=scale * jn,
            I1=scale * (x1 * jn + d * prod(y)),
            I2=scale * (x1^2 * jn + 2 * x1 * d * prod(y) + d^2 * second(y)),
        )
    end
    function pole(h, r)
        y = y_of_r(r)
        p = (x1 - h) / (x2 - h)
        alpha = (1 - n) / (p - n)
        beta = (p - 1) / (p - n)
        return scale / (x2 - h) * (alpha * _j_inv(n, y) + beta * _j_inv(p, y))
    end
    function inverse_i0(target)
        jn = max(target / scale, 0.0)
        z = min(tanh(sqrt(n) * jn), prevfloat(1.0))
        y = z / sqrt(n)
        return (x2 - x1 * y^2) / (1 - y^2)
    end
    return (kind=:k11_hyperbolic_outer_double, roots=(x1=x1, x2=x2, x3=x3),
        upper=x3, basis=basis, pole=pole,
        inverse_i0=inverse_i0)
end

function _k7_model(roots)
    x1, x3 = roots
    d = x3 - x1
    scale = sqrt(2 / d)
    u_of_r(r) = sqrt(max((r - x1) / d, 1.0))
    j0(u) = 0.5 * log((u - 1) / (u + 1))
    function basis(r)
        u = u_of_r(r)
        return (
            I0=scale * j0(u),
            I1=scale * (x3 * j0(u) + d * u),
            I2=scale * (x3^2 * j0(u) + 2 * x3 * d * u +
                d^2 * (u^3 / 3 - u)),
        )
    end
    function pole(h, r)
        u = u_of_r(r)
        return scale * (j0(u) / (x3 - h) - d / (x3 - h) *
            _quadratic_denominator_primitive(x1 - h, d, u))
    end
    function inverse_i0(target)
        target < 0 || return Inf
        u = -inv(tanh(target / scale))
        return x1 + d * u^2
    end
    return (
        kind=:k7_parabolic_exterior_repeated_root,
        roots=(x1=x1, x3=x3),
        lower=x3,
        basis=basis,
        pole=pole,
        inverse_i0=inverse_i0,
        infinity_i0=0.0,
    )
end

function _k10_model(energy, roots)
    x1, x2, x3 = roots
    lead = _e2m1(energy)
    k = (x3 - x2) / (x3 - x1)
    d = x2 - x1
    h0 = d * k
    scale = 2 / sqrt(lead * (x3 - x1) * (x3 - x2))
    z_of_r(r) = sqrt(clamp(
        (r - x1) * (x3 - x2) / ((r - x2) * (x3 - x1)),
        k,
        1.0,
    ))
    j0(z) = atanh(z)
    jk(z) = _j_z2_minus(k, z)
    jk2(z) = -z / (2 * k * (z^2 - k)) - jk(z) / (2 * k)
    prod(z) = (j0(z) + jk(z)) / (1 - k)
    second(z) = (j0(z) + jk(z)) / (1 - k)^2 + jk2(z) / (1 - k)
    function basis(r)
        z = z_of_r(r)
        return (
            I0=-scale * j0(z),
            I1=-scale * (x2 * j0(z) + h0 * prod(z)),
            I2=-scale * (x2^2 * j0(z) + 2 * x2 * h0 * prod(z) +
                h0^2 * second(z)),
        )
    end
    function pole(h, r)
        z = z_of_r(r)
        kh = k * (x1 - h) / (x2 - h)
        alpha = (1 - k) / (1 - kh)
        beta = (kh - k) / (1 - kh)
        return -scale / (x2 - h) * (alpha * j0(z) +
            beta * _j_z2_minus(kh, z))
    end
    function inverse_i0(target)
        z = tanh(max(-target / scale, 0.0))
        z = clamp(z, nextfloat(sqrt(k)), prevfloat(1.0))
        return x2 + h0 / (z^2 - k)
    end
    return (
        kind=:k10_hyperbolic_exterior_repeated_root,
        roots=(x1=x1, x2=x2, x3=x3),
        lower=x3,
        basis=basis,
        pole=pole,
        inverse_i0=inverse_i0,
        infinity_i0=-scale * atanh(sqrt(k)),
    )
end

function _homoclinic_model(energy, x1, rc, ra)
    kappa = -_e2m1(energy)
    kappa > 0 || error("The homoclinic model requires E<1.")
    x1 < rc < ra || error("Homoclinic roots must satisfy x1<rc<ra.")
    total = ra - x1
    outer_gap = ra - rc
    alpha = sqrt((rc - x1) / outer_gap)
    scale = 2 / (sqrt(kappa) * outer_gap)

    function t_of_r(r)
        rc < r <= ra || throw(DomainError(
            r, "A homoclinic radius must lie in (rc,ra]."))
        return sqrt((r - x1) / (ra - r))
    end
    function t_of_inner_r(r)
        x1 < r < rc || throw(DomainError(
            r, "A whirling radius must lie in (x1,rc)."))
        return sqrt((r - x1) / (ra - r))
    end
    jalpha(t) = log(abs((t - alpha) / (t + alpha))) / (2 * alpha)
    j1(t) = atan(t) - pi / 2
    j2(t) = 0.5 * (atan(t) - pi / 2 + t / (1 + t^2))
    a1(t) = (jalpha(t) - j1(t)) / (1 + alpha^2)
    a2(t) = (jalpha(t) - j1(t)) / (1 + alpha^2)^2 -
        j2(t) / (1 + alpha^2)

    function basis(r)
        if r == ra
            return (I0=0.0, I1=0.0, I2=0.0)
        end
        t = t_of_r(r)
        ja = jalpha(t)
        first = a1(t)
        second = a2(t)
        return (
            I0=scale * ja,
            I1=scale * (ra * ja - total * first),
            I2=scale * (ra^2 * ja - 2 * ra * total * first +
                total^2 * second),
        )
    end

    function inner_basis(r)
        t = t_of_inner_r(r)
        ja = jalpha(t)
        first = a1(t)
        second = a2(t)
        return (
            I0=-scale * ja,
            I1=-scale * (ra * ja - total * first),
            I2=-scale * (ra^2 * ja - 2 * ra * total * first +
                total^2 * second),
        )
    end

    function pole(h, r)
        r == ra && return 0.0
        t = t_of_r(r)
        a0 = x1 - h
        b0 = ra - h
        coefficient_repeated = 1 / (rc - h)
        coefficient_quadratic = (rc - ra) / (rc - h)
        quadratic_infinity = a0 > 0 ? pi / (2 * sqrt(a0 * b0)) : 0.0
        return scale * (
            coefficient_repeated * jalpha(t) +
            coefficient_quadratic *
                (_quadratic_denominator_primitive(a0, b0, t) -
                 quadratic_infinity)
        )
    end

    function inner_pole(h, r)
        t = t_of_inner_r(r)
        a0 = x1 - h
        b0 = ra - h
        coefficient_repeated = 1 / (rc - h)
        coefficient_quadratic = (rc - ra) / (rc - h)
        quadratic_infinity = a0 > 0 ? pi / (2 * sqrt(a0 * b0)) : 0.0
        return -scale * (
            coefficient_repeated * jalpha(t) +
            coefficient_quadratic *
                (_quadratic_denominator_primitive(a0, b0, t) -
                 quadratic_infinity)
        )
    end

    function inverse_i0(target)
        target <= 0 || throw(DomainError(
            target, "The homoclinic I0 target must not exceed zero."))
        target == 0 && return ra
        exponential = exp(2 * alpha * target / scale)
        t = alpha * (1 + exponential) / (1 - exponential)
        return (x1 + ra * t^2) / (1 + t^2)
    end


    function inner_inverse_i0(target)
        target >= 0 || throw(DomainError(
            target, "The whirling I0 target must be nonnegative."))
        t = alpha * tanh(alpha * target / scale)
        return (x1 + ra * t^2) / (1 + t^2)
    end

    return (
        kind=:homoclinic_outer_simple_inner_double,
        roots=(x1=x1, rc=rc, ra=ra),
        lower=rc,
        upper=ra,
        basis=basis,
        pole=pole,
        inverse_i0=inverse_i0,
        inner_basis=inner_basis,
        inner_pole=inner_pole,
        inner_inverse_i0=inner_inverse_i0,
        scale=scale,
        alpha=alpha,
    )
end

function _j_z2_minus(k, z)
    k > 0 || error("The quadratic parameter must be positive.")
    root = sqrt(k)
    if z < root
        return -atanh(z / root) / root
    elseif z > root
        return log(abs((z - root) / (z + root))) / (2 * root)
    end
    return -Inf
end
