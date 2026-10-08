# Radial models for repeated roots inside the horizon (B7-B9, C6-C12).

function _qf(a, b, y)
    iszero(b) && return y / a
    iszero(a) && return -1 / (b * y)
    if b > 0
        return _quadratic_denominator_primitive(a, b, y)
    elseif a > 0
        return atanh(y * sqrt(-b / a)) / sqrt(-a * b)
    end
    error("The elementary quadratic primitive is outside its physical branch.")
end

# _qf(a, b, y) − _qf(a, b, 1) given 1 − y (branches as _qf), without subtracting the two
# primitives: atan and atanh differences by their addition formulas, the log one as a log1p
function _qf_from_one(a, b, y, omy)
    iszero(b) && return -omy / a
    iszero(a) && return -omy / (b * y)
    if a > 0 && b > 0
        t = sqrt(b / a)
        return -atan(omy * t / (1 + y * t^2)) / sqrt(a * b)
    elseif a > 0 && b < 0
        t = sqrt(-b / a)
        return -atanh(omy * t / (1 - y * t^2)) / sqrt(-a * b)
    elseif a < 0 && b > 0
        c = sqrt(-a / b)
        return log1p(-2c * omy / ((y + c) * (1 - c))) / (2 * sqrt(-a * b))
    end
    error("The elementary quadratic primitive is outside its physical branch.")
end

function _qfinv(a, b, value)
    iszero(b) && return a * value
    iszero(a) && return -1 / (b * value)
    if a > 0 && b > 0
        product = a * b
        scale = sqrt(product)
        scale_error = (fma(a, b, -product) + fma(-scale, scale, product)) / (2scale)
        phase = scale * value
        phase_error = fma(scale, value, -phase) + scale_error * value
        # Preserve the phase residual before the tangent amplifies it near a pole.
        tangent = tan(phase)
        correction = tan(phase_error)
        return sqrt(a / b) * (tangent + correction) / (1 - tangent * correction)
    elseif a > 0 && b < 0
        return sqrt(a / -b) * tanh(sqrt(-a * b) * value)
    elseif a < 0 && b > 0
        return -sqrt(-a / b) / tanh(sqrt(-a * b) * value)
    end
    error("The elementary quadratic inverse is outside its physical branch.")
end

function _hf(c, eta, w)
    d = hypot(c, eta)
    return log(abs((c * w + eta - d) /
        (c * w + eta + d))) / (2d)
end

function _interior_double_model(case_id, a, energy, lz, d, s, outer)
    kappa = -_e2m1(energy)
    aa = outer - d
    bb = s - d
    y(r) = sqrt(max((outer - r) / (r - s), 0.0))
    radius(yvalue) = (outer + s * yvalue^2) / (1 + yvalue^2)
    function raw(r)
        yy = y(r)
        f0 = _qf(aa, bb, yy)
        i0 = 2 * f0 / sqrt(kappa)
        i1 = d * i0 + 2atan(yy) / sqrt(kappa)
        i2 = d^2 * i0 + (4d + aa + bb) * atan(yy) / sqrt(kappa) +
            (aa - bb) * yy / (sqrt(kappa) * (1 + yy^2))
        return (I0=i0, I1=i1, I2=i2)
    end
    basis(r) = begin
        value = raw(r)
        (I0=-value.I0, I1=-value.I1, I2=-value.I2)
    end
    pole(h, r) = begin
        yy = y(r)
        -2 * (_qf(aa, bb, yy) - _qf(outer - h, s - h, yy)) /
            (sqrt(kappa) * (d - h))
    end
    i0_turn = basis(outer).I0
    radius_at(δ) = radius(_qfinv(aa, bb, -(i0_turn - δ) * sqrt(kappa) / 2))
    mino(r) = i0_turn - basis(r).I0
    return (
        kind=case_id === :B7 ? :elliptic_double_d_below_simple :
            :elliptic_double_simple_below_d,
        roots=(repeated=d, simple=s, outer=outer), turn=outer,
        basis=basis, pole=pole, radius=radius_at, mino=mino, inward=true,
    )
end

function _parabolic_double_model(case_id, a, energy, lz, d, s)
    rt2 = sqrt(_float_type(a, energy, lz, d, s)(2))
    bb = s - d
    y(r) = sqrt(max(r - s, 0.0))
    # I0 is measured from infinity: for bb > 0, atan(y/√bb) − π/2 = −atan(√bb/y)
    function basis(r)
        yy = y(r)
        i0 = bb > 0 ? -sqrt(2 / bb) * atan(sqrt(bb) / yy) : rt2 * _qf(bb, 1.0, yy)
        i1 = d * i0 + rt2 * yy
        i2 = d^2 * i0 + rt2 *
            ((2d + bb) * yy + yy^3 / 3)
        return (I0=i0, I1=i1, I2=i2)
    end
    pole(h, r) = begin
        yy = y(r)
        rt2 * (_qf(bb, 1.0, yy) - _qf(s - h, 1.0, yy)) /
            (d - h)
    end
    # r at Mino time δ after the infinity endpoint: y = √bb cot(√bb δ/√2) for bb > 0,
    # √(−bb) coth(√(−bb) δ/√2) for bb < 0
    radius(δ) = s + (bb > 0 ? bb / tan(sqrt(bb) * δ / rt2)^2 :
        -bb / tanh(sqrt(-bb) * δ / rt2)^2)
    return (
        kind=case_id === :C6 ? :parabolic_double_d_below_simple :
            :parabolic_double_simple_below_d,
        roots=(repeated=d, simple=s), lower=max(d, s),
        basis=basis, pole=pole, radius=radius, mino=r -> -basis(r).I0,
        inward=true,
    )
end

function _hyperbolic_real_double_model(case_id, a, energy, lz,
        d, s1, s2)
    lead = _e2m1(energy)
    aa = s2 - d
    bb = d - s1
    y(r) = sqrt(max((r - s2) / (r - s1), 0.0))
    # I0 is measured from infinity (y = 1): 1 − y = (s2 − s1)/((r − s1)(1 + y))
    function basis(r)
        yy = y(r)
        i0 = 2 * _qf_from_one(aa, bb, yy, (s2 - s1) / ((r - s1) * (1 + yy))) / sqrt(lead)
        at = atanh(yy)
        i1 = d * i0 + 2at / sqrt(lead)
        i2 = d^2 * i0 + (4d + aa - bb) * at / sqrt(lead) +
            (aa + bb) * yy / (sqrt(lead) * (1 - yy^2))
        return (I0=i0, I1=i1, I2=i2)
    end
    pole(h, r) = begin
        yy = y(r)
        2 * (_qf(aa, bb, yy) - _qf(s2 - h, h - s1, yy)) /
            (sqrt(lead) * (d - h))
    end
    # r = (s2 − s1 y²)/(1 − y²) at Mino time δ after the infinity endpoint, where y → 1. In the
    # angle θ of _qfinv (y = √(a/b) tan θ, √(a/−b) tanh θ or √(−a/b) coth θ, y = 1 at θ∞),
    # 1 − y² is a product of the offset η = √|ab| √lead δ/2 and θ∞ ± θ; y² itself is taken
    # from θ, so that neither the far field nor the horizon side cancels
    function radius(δ)
        η = sqrt(abs(aa * bb) * lead) * δ / 2
        if aa > 0 && bb > 0
            θ = atan(sqrt(bb / aa))
            y2 = aa / bb * tan(θ - η)^2
            gap = (aa + bb) / bb * sin(η) * sin(2θ - η) / cos(θ - η)^2
        elseif aa > 0
            θ = atanh(sqrt(-bb / aa))
            y2 = aa / -bb * tanh(θ - η)^2
            gap = (aa + bb) / -bb * sinh(η) * sinh(2θ - η) / cosh(θ - η)^2
        else
            θ = atanh(sqrt(-aa / bb))
            y2 = -aa / bb / tanh(θ + η)^2
            gap = (aa + bb) / bb * sinh(η) * sinh(2θ + η) / sinh(θ + η)^2
        end
        return (s2 - s1 * y2) / gap
    end
    return (
        kind=case_id === :C9 ? :hyperbolic_simple_double_simple :
            :hyperbolic_two_simple_below_double,
        roots=(simple1=s1, repeated=d, simple2=s2), lower=max(d, s2),
        basis=basis, pole=pole, radius=radius, mino=r -> -basis(r).I0,
        inward=true,
    )
end

function _hyperbolic_complex_double_model(a, energy, lz, d, rho, eta)
    lead = _e2m1(energy)
    c = d - rho
    D = hypot(c, eta)
    σ = c + eta + D
    w(r) = (r - rho) / (hypot(r - rho, eta) + eta)
    # 1 − w without cancellation: h − (r − ρ) = η²/(h + r − ρ) for r > ρ, h = √((r − ρ)² + η²)
    function one_minus_w(r)
        h = hypot(r - rho, eta)
        return (eta + (r > rho ? eta^2 / (h + r - rho) : h - (r - rho))) / (h + eta)
    end
    # I0 is measured from infinity (w = 1): _hf(w) − _hf(1) = log(g/g∞)/(2D) with
    # g/g∞ − 1 = (1 − w)(2cη − σ²)/(2η(cw + η + D)), σ = c + η + D
    function basis(r)
        ww = w(r)
        u = 2atanh(ww)
        i0 = log1p(one_minus_w(r) * (2c * eta - σ^2) / (2eta * (c * ww + eta + D))) /
            (D * sqrt(lead))
        i1 = d * i0 + u / sqrt(lead)
        i2 = d^2 * i0 + (2d - c) * u / sqrt(lead) +
            eta * cosh(u) / sqrt(lead)
        return (I0=i0, I1=i1, I2=i2)
    end
    pole(h, r) = begin
        ww = w(r)
        2 * (_hf(c, eta, ww) - _hf(h - rho, eta, ww)) /
            (sqrt(lead) * (d - h))
    end
    # r = ρ + 2ηw/(1 − w²) at Mino time δ after the infinity endpoint (w → 1). With
    # g = (cw + η − D)/(cw + η + D), D = √(c² + η²), the primitive gives g = g∞ e^(−x),
    # x = D√lead δ, and g∞ = 2cη/(c + η + D)²; then 1 − w = (z∞ − z)/c with z = D(1 + g)/(1 − g)
    # is 4Dη(1 − e^(−x))/((c + η + D)²(1 − g∞)(1 − g)), free of cancellation
    σ2 = σ^2
    g∞ = 2c * eta / σ2
    function radius(δ)
        x = D * sqrt(lead) * δ
        g = g∞ * exp(-x)
        gap = 4D * eta * -expm1(-x) / (σ2 * (1 - g∞) * (1 - g))
        return rho + 2eta * (1 - gap) / (gap * (2 - gap))
    end
    return (
        kind=:hyperbolic_double_plus_complex_pair,
        roots=(repeated=d, rho=rho, eta=eta), lower=d,
        basis=basis, pole=pole, radius=radius, mino=r -> -basis(r).I0,
        inward=true,
    )
end

function _interior_triple_model(a, energy, lz, d, outer)
    kappa = -_e2m1(energy)
    c = 2 / (sqrt(kappa) * (outer - d))
    y(r) = sqrt(max((outer - r) / (r - d), 0.0))
    radius(yvalue) = (outer + d * yvalue^2) / (1 + yvalue^2)
    function raw(r)
        yy = y(r)
        i0 = c * yy
        i1 = c * (d * yy + (outer - d) * atan(yy))
        i2 = c * (d^2 * yy + 2d * (outer - d) * atan(yy) +
            (outer - d)^2 * (atan(yy) + yy / (1 + yy^2)) / 2)
        return (I0=i0, I1=i1, I2=i2)
    end
    basis(r) = begin
        value = raw(r)
        (I0=-value.I0, I1=-value.I1, I2=-value.I2)
    end
    pole(h, r) = begin
        yy = y(r)
        -c * (yy / (d - h) +
            (d - outer) * _qf(outer - h, d - h, yy) / (d - h))
    end
    i0_turn = basis(outer).I0
    radius_at(δ) = radius(-(i0_turn - δ) / c)
    mino(r) = i0_turn - basis(r).I0
    return (
        kind=:elliptic_triple_below_horizon,
        roots=(repeated=d, outer=outer), turn=outer,
        basis=basis, pole=pole, radius=radius_at, mino=mino, inward=true,
    )
end

function _parabolic_triple_model(a, energy, lz, d)
    rt2 = sqrt(_float_type(a, energy, lz, d)(2))
    v(r) = inv(sqrt(r - d))
    function basis(r)
        vv = v(r)
        return (
            I0=-rt2 * vv,
            I1=-rt2 * (d * vv - inv(vv)),
            I2=-rt2 * (d^2 * vv - 2d / vv - inv(3vv^3)),
        )
    end
    pole(h, r) = begin
        vv = v(r)
        -rt2 * (vv - _qf(1.0, d - h, vv)) / (d - h)
    end
    # r at Mino time δ after the infinity endpoint (where I0 = 0): v = δ/√2
    radius(δ) = d + 2 / δ^2
    return (
        kind=:parabolic_triple_below_horizon,
        roots=(repeated=d,), lower=d,
        basis=basis, pole=pole, radius=radius, mino=r -> -basis(r).I0,
        inward=true,
    )
end

function _hyperbolic_triple_model(a, energy, lz, s, d)
    lead = _e2m1(energy)
    b = d - s
    c = 2 / (sqrt(lead) * b)
    v(r) = sqrt((r - s) / (r - d))
    function basis(r)
        vv = v(r)
        g = log((vv - 1) / (vv + 1)) / 2
        h2 = -vv / (2 * (vv^2 - 1)) - g / 2
        return (
            I0=-c * b / ((r - d) * (vv + 1)),        # −c(v − 1), measured from infinity
            I1=-c * (d * vv + b * g),
            I2=-c * (d^2 * vv + 2d * b * g + b^2 * h2),
        )
    end
    pole(h, r) = begin
        vv = v(r)
        -c * (vv / (d - h) -
            b * _qf(h - s, d - h, vv) / (d - h))
    end
    # r = d + (d − s)/(v² − 1) at Mino time δ after the infinity endpoint: v − 1 = δ/c exactly
    radius(δ) = (e = δ / c; d + (d - s) / (e * (2 + e)))
    return (
        kind=:hyperbolic_triple_below_horizon,
        roots=(simple=s, repeated=d), lower=d,
        basis=basis, pole=pole, radius=radius, mino=r -> -basis(r).I0,
        inward=true,
    )
end

function _root_data(structure)
    items = collect(structure.real_roots)
    return items, _nonreal_roots(structure)
end

"""Return the elementary radial model for an admitted interior repeated root."""
function interior_repeated_radial_model(case_id, a, energy, lz, q, structure)
    items, nonreal = _root_data(structure)
    repeated = only(filter(item -> item.multiplicity >= 2, items)).radius
    simple = sort([item.radius for item in items if item.multiplicity == 1])
    if case_id in (:B7, :B8)
        return _interior_double_model(
            case_id, a, energy, lz, repeated, simple[1], simple[2])
    elseif case_id === :B9
        return _interior_triple_model(
            a, energy, lz, repeated, only(simple))
    elseif case_id in (:C6, :C7)
        return _parabolic_double_model(
            case_id, a, energy, lz, repeated, only(simple))
    elseif case_id === :C8
        return _parabolic_triple_model(a, energy, lz, repeated)
    elseif case_id in (:C9, :C10)
        return _hyperbolic_real_double_model(
            case_id, a, energy, lz, repeated, simple[1], simple[2])
    elseif case_id === :C11
        length(nonreal) == 2 || error("C11 requires one conjugate root pair.")
        upper = nonreal[argmax(imag.(nonreal))]
        return _hyperbolic_complex_double_model(
            a, energy, lz, repeated, real(upper), abs(imag(upper)))
    elseif case_id === :C12
        return _hyperbolic_triple_model(
            a, energy, lz, only(simple), repeated)
    end
    error("No interior repeated-root radial model for $(case_id).")
end
