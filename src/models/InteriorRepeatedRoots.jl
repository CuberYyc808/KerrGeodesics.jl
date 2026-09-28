# Radial models for repeated roots inside the horizon (B7-B9, C6-C12).

function _qf(a, b, y)
    scale = max(1.0, abs(a), abs(b))
    abs(b) <= 32eps(Float64) * scale && return y / a
    abs(a) <= 32eps(Float64) * scale && return -1 / (b * y)
    if a > 0 && b > 0
        return atan(y * sqrt(b / a)) / sqrt(a * b)
    elseif a > 0 && b < 0
        return atanh(y * sqrt(-b / a)) / sqrt(-a * b)
    elseif a < 0 && b > 0
        c = sqrt(-a / b)
        return log(abs((y - c) / (y + c))) / (2 * sqrt(-a * b))
    end
    error("The elementary quadratic primitive is outside its physical branch.")
end

function _qfinv(a, b, value)
    scale = max(1.0, abs(a), abs(b))
    abs(b) <= 32eps(Float64) * scale && return a * value
    abs(a) <= 32eps(Float64) * scale && return -1 / (b * value)
    if a > 0 && b > 0
        return sqrt(a / b) * tan(sqrt(a * b) * value)
    elseif a > 0 && b < 0
        return sqrt(a / -b) * tanh(sqrt(-a * b) * value)
    elseif a < 0 && b > 0
        exponent = exp(2 * sqrt(-a * b) * value)
        return sqrt(-a / b) * (1 + exponent) / (1 - exponent)
    end
    error("The elementary quadratic inverse is outside its physical branch.")
end

function _hf(c, eta, w)
    d = hypot(c, eta)
    return log(abs((c * w + eta - d) /
        (c * w + eta + d))) / (2d)
end

function _hfinv(c, eta, value)
    d = hypot(c, eta)
    exponent = exp(2d * value)
    zinner = d * (1 - exponent) / (1 + exponent)
    zouter = d * (1 + exponent) / (1 - exponent)
    candidates = ((zinner - eta) / c, (zouter - eta) / c)
    matches = filter(w -> abs(w) < 1 + 1e-12 &&
        isapprox(_hf(c, eta, w), value; atol=2e-11, rtol=2e-11),
        candidates)
    length(matches) == 1 || error(
        "The elementary hyperbolic-log inverse branch is ambiguous.")
    return only(matches)
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
    inverse_i0(target) = radius(_qfinv(
        aa, bb, -target * sqrt(kappa) / 2))
    return (
        kind=case_id === :B7 ? :elliptic_double_d_below_simple :
            :elliptic_double_simple_below_d,
        roots=(repeated=d, simple=s, outer=outer), upper=outer,
        basis=basis, pole=pole,
        inverse_i0=inverse_i0,
    )
end

function _parabolic_double_model(case_id, a, energy, lz, d, s)
    bb = s - d
    y(r) = sqrt(max(r - s, 0.0))
    radius(yvalue) = s + yvalue^2
    function basis(r)
        yy = y(r)
        f0 = _qf(bb, 1.0, yy)
        i0 = sqrt(2.0) * f0
        i1 = d * i0 + sqrt(2.0) * yy
        i2 = d^2 * i0 + sqrt(2.0) *
            ((2d + bb) * yy + yy^3 / 3)
        return (I0=i0, I1=i1, I2=i2)
    end
    pole(h, r) = begin
        yy = y(r)
        sqrt(2.0) * (_qf(bb, 1.0, yy) - _qf(s - h, 1.0, yy)) /
            (d - h)
    end
    inverse_i0(target) = radius(_qfinv(bb, 1.0, target / sqrt(2.0)))
    infinity_i0 = bb > 0 ? pi / sqrt(2bb) : 0.0
    return (
        kind=case_id === :C6 ? :parabolic_double_d_below_simple :
            :parabolic_double_simple_below_d,
        roots=(repeated=d, simple=s), lower=max(d, s),
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=infinity_i0,
    )
end

function _hyperbolic_real_double_model(case_id, a, energy, lz,
        d, s1, s2)
    lead = _e2m1(energy)
    aa = s2 - d
    bb = d - s1
    y(r) = sqrt(max((r - s2) / (r - s1), 0.0))
    radius(yvalue) = (s2 - s1 * yvalue^2) / (1 - yvalue^2)
    function basis(r)
        yy = y(r)
        f0 = _qf(aa, bb, yy)
        i0 = 2 * f0 / sqrt(lead)
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
    inverse_i0(target) = radius(_qfinv(
        aa, bb, target * sqrt(lead) / 2))
    infinity_i0 = 2 * _qf(aa, bb, 1.0) / sqrt(lead)
    return (
        kind=case_id === :C9 ? :hyperbolic_simple_double_simple :
            :hyperbolic_two_simple_below_double,
        roots=(simple1=s1, repeated=d, simple2=s2), lower=max(d, s2),
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=infinity_i0,
    )
end

function _hyperbolic_complex_double_model(a, energy, lz, d, rho, eta)
    lead = _e2m1(energy)
    c = d - rho
    w(r) = (r - rho) / (hypot(r - rho, eta) + eta)
    radius(wvalue) = rho + 2eta * wvalue / (1 - wvalue^2)
    function basis(r)
        ww = w(r)
        u = 2atanh(ww)
        f0 = _hf(c, eta, ww)
        i0 = 2 * f0 / sqrt(lead)
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
    inverse_i0(target) = radius(_hfinv(
        c, eta, target * sqrt(lead) / 2))
    infinity_i0 = 2 * _hf(c, eta, 1.0) / sqrt(lead)
    return (
        kind=:hyperbolic_double_plus_complex_pair,
        roots=(repeated=d, rho=rho, eta=eta), lower=d,
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=infinity_i0,
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
    inverse_i0(target) = radius(-target / c)
    return (
        kind=:elliptic_triple_below_horizon,
        roots=(repeated=d, outer=outer), upper=outer,
        basis=basis, pole=pole,
        inverse_i0=inverse_i0,
    )
end

function _parabolic_triple_model(a, energy, lz, d)
    v(r) = inv(sqrt(r - d))
    radius(vvalue) = d + inv(vvalue^2)
    function basis(r)
        vv = v(r)
        return (
            I0=-sqrt(2.0) * vv,
            I1=-sqrt(2.0) * (d * vv - inv(vv)),
            I2=-sqrt(2.0) * (d^2 * vv - 2d / vv - inv(3vv^3)),
        )
    end
    pole(h, r) = begin
        vv = v(r)
        -sqrt(2.0) * (vv - _qf(1.0, d - h, vv)) / (d - h)
    end
    inverse_i0(target) = radius(-target / sqrt(2.0))
    return (
        kind=:parabolic_triple_below_horizon,
        roots=(repeated=d,), lower=d,
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=0.0,
    )
end

function _hyperbolic_triple_model(a, energy, lz, s, d)
    lead = _e2m1(energy)
    b = d - s
    c = 2 / (sqrt(lead) * b)
    v(r) = sqrt((r - s) / (r - d))
    radius(vvalue) = (d * vvalue^2 - s) / (vvalue^2 - 1)
    function basis(r)
        vv = v(r)
        g = log((vv - 1) / (vv + 1)) / 2
        h2 = -vv / (2 * (vv^2 - 1)) - g / 2
        return (
            I0=-c * vv,
            I1=-c * (d * vv + b * g),
            I2=-c * (d^2 * vv + 2d * b * g + b^2 * h2),
        )
    end
    pole(h, r) = begin
        vv = v(r)
        -c * (vv / (d - h) -
            b * _qf(h - s, d - h, vv) / (d - h))
    end
    inverse_i0(target) = radius(-target / c)
    return (
        kind=:hyperbolic_triple_below_horizon,
        roots=(simple=s, repeated=d), lower=d,
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=-c,
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
