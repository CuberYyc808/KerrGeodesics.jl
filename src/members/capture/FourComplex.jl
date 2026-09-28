# The capture C5 (E > 1, four complex radial roots): the real Jacobi form of r(λ) for two
# complex-conjugate root pairs, its Mino-time basis and horizon poles, and the constructor.

# the Jacobi parameters of R = (E² − 1)|r − u − iw|²|r − v − is|² (classified constants)
function _four_complex_parameters(energy, structure)
    upper_roots = sort([root for root in structure.raw_roots if imag(root) > 0];
        by=real)
    length(upper_roots) == 2 || error(
        "Could not resolve two upper-half-plane radial roots.")
    u = real(upper_roots[1])
    w = imag(upper_roots[1])
    v = real(upper_roots[2])
    s = imag(upper_roots[2])
    delta = u - v
    scale = max(1.0, abs(u), abs(v), w, s)
    delta < -128 * eps(Float64) * scale || error(
        "The C5 Jacobi form requires distinct real parts of the two complex radial-root pairs.")
    w > 0 && s > 0 || error("Complex radial-root heights must be positive.")

    d0 = delta^2 + w^2 + s^2
    discriminant = d0^2 - 4 * s^2 * w^2
    discriminant >= -256 * eps(Float64) * max(1.0, d0^2) || error(
        "Four-complex-root Jacobi discriminant is negative.")
    root = sqrt(max(discriminant, 0.0))
    lambda1 = (d0 + root) / (2 * s^2)
    lambda2 = (d0 - root) / (2 * s^2)
    lambda1 > 1 > lambda2 > 0 || error(
        "Four-complex-root Jacobi parameters require lambda1>1>lambda2>0.")
    c = lambda1 - 1
    d = lambda1 - lambda2
    c < d < c + 1 || error("Four-complex-root Jacobi ordering failed.")
    m = d / lambda1
    0 < m < 1 || error("Four-complex-root Jacobi modulus must lie in (0,1).")
    lead = _e2m1(energy)
    alpha = v - delta / c
    omega = s * sqrt(lambda1 * lead)
    n = d / c
    zinf = sqrt(c / d)

    relation_residual = maximum(abs.((
        delta^2 - s^2 * c * (d - c),
        w^2 - s^2 * lambda1 * lambda2,
    ))) / max(1.0, delta^2, w^2, s^2)
    relation_residual <= 5.0e-11 || error(
        "The four-complex-root Jacobi parameters violate their defining identities (relative residual $(relation_residual) > 5e-11).")
    return (
        u=u, w=w, v=v, s=s, delta=delta,
        lambda1=lambda1, lambda2=lambda2, c=c, d=d, m=m,
        lead=lead, alpha=alpha, omega=omega, n=n, zinf=zinf,
        relation_residual=relation_residual,
    )
end

function _four_complex_partial_fractions(numerator, denominator)
    n0, n1coef, n2coef = numerator
    p0, p1, p2 = denominator
    scale = max(1.0, abs(p0), abs(p1), abs(p2))
    abs(p0) > 128 * eps(Float64) * scale || error(
        "Four-complex-root pole transform is at its zero-constant limit.")
    p2 > 0 || error("Four-complex-root pole quadratic must have positive leading coefficient.")
    discriminant = p1^2 - 4 * p2 * p0
    tolerance = 512 * eps(Float64) * scale^2
    discriminant >= -tolerance || error(
        "Four-complex-root pole characteristics are not real.")
    root = sqrt(max(discriminant, 0.0))
    quotient = n2coef / p2
    remainder0 = n0 - quotient * p0
    remainder1 = n1coef - quotient * p1

    if root <= sqrt(tolerance)
        characteristic = -p1 / (2 * p0)
        abs(characteristic) > 128 * eps(Float64) || error(
            "Repeated pole characteristic is zero.")
        first = -remainder1 / (p0 * characteristic)
        second = remainder0 / p0 - first
        return (kind=:repeated, constant=quotient,
            characteristics=(characteristic,),
            coefficients=(first, second), discriminant=discriminant)
    end

    x1 = (-p1 - root) / (2 * p2)
    x2 = (-p1 + root) / (2 * p2)
    x1 > 0 && x2 > 0 || error(
        "Four-complex-root pole roots must be positive in z squared.")
    characteristic1 = inv(x1)
    characteristic2 = inv(x2)
    abs(characteristic1 - characteristic2) > 128 * eps(Float64) *
        max(1.0, abs(characteristic1), abs(characteristic2)) || error(
        "Distinct pole characteristics collapsed numerically.")
    sum0 = remainder0 / p0
    sum1 = remainder1 / p0
    first = (sum1 + characteristic1 * sum0) /
        (characteristic1 - characteristic2)
    second = sum0 - first
    return (kind=:distinct, constant=quotient,
        characteristics=(characteristic1, characteristic2),
        coefficients=(first, second), discriminant=discriminant)
end

function _four_complex_even_primitive(data, z, m)
    phi = asin(clamp(z, -1.0, 1.0))
    value = data.constant * Elliptic.F(phi, m)
    if data.kind === :distinct
        for (characteristic, coefficient) in
                zip(data.characteristics, data.coefficients)
            value += coefficient * _pi_real(characteristic, phi, m)
        end
    else
        characteristic = only(data.characteristics)
        value += data.coefficients[1] * _pi_real(characteristic, phi, m)
        value += data.coefficients[2] *
            _j2_legendre(characteristic, m, phi)
    end
    return value
end

function _four_complex_odd_primitive(data, z, m)
    abs(data.constant) <= 64 * eps(Float64) || error(
        "Odd four-complex-root pole decomposition acquired a polynomial part.")
    if data.kind === :distinct
        return sum(coefficient * _odd_h1(
            characteristic, m, z) for
            (characteristic, coefficient) in
                zip(data.characteristics, data.coefficients))
    end
    characteristic = only(data.characteristics)
    return data.coefficients[1] * _odd_h1(
        characteristic, m, z) +
        data.coefficients[2] * _odd_h2(
            characteristic, m, z)
end

function _four_complex_model(parameters)
    p = parameters
    z_of_r(r) = isinf(r) ? p.zinf :
        sqrt(p.c / p.d) * (r - p.alpha) / hypot(r - p.v, p.s)

    function radius_from_z(z)
        zv = float(z)
        abs(zv - p.zinf) <= 32 * eps(Float64) && return Inf
        abs(zv) < 1 || throw(DomainError(z, "Four-complex-root Jacobi z must lie in (-1,1)."))
        denominator = p.d * zv^2 - p.c
        abs(denominator) > 64 * eps(Float64) * max(1.0, p.c) || return Inf
        numerator = p.delta - p.s * p.d * zv * sqrt(max(1 - zv^2, 0.0))
        return p.v + numerator / denominator
    end

    function basis(r)
        z = z_of_r(r)
        phi = asin(clamp(z, -1.0, 1.0))
        f = Elliptic.F(phi, p.m)
        pin = _pi_real(p.n, phi, p.m)
        j2 = _j2_legendre(p.n, p.m, phi)
        h1 = _odd_h1(p.n, p.m, z)
        h2 = _odd_h2(p.n, p.m, z)
        acoef = -p.delta / p.c
        bcoef = p.s * p.d / p.c
        rat_a = -inv(p.n^2)
        rat_b = 2 / p.n^2 - inv(p.n)
        rat_c = (p.n - 1) / p.n^2
        return (
            I0=f / p.omega,
            I1=(p.v * f + acoef * pin + bcoef * h1) / p.omega,
            I2=(p.v^2 * f + 2 * p.v * acoef * pin +
                2 * p.v * bcoef * h1 + acoef^2 * j2 +
                2 * acoef * bcoef * h2 + bcoef^2 *
                (rat_a * f + rat_b * pin + rat_c * j2)) / p.omega,
        )
    end

    function pole_data(h)
        hgap = p.v - h
            aa = hgap * p.d
            bb = p.delta - hgap * p.c
            ss = p.s * p.d
            denominator = (
                bb^2,
                2 * aa * bb - ss^2,
                aa^2 + ss^2,
            )
            even = _four_complex_partial_fractions((
                -p.c * bb,
                p.d * bb - p.c * aa,
                p.d * aa,
            ), denominator)
            odd = _four_complex_partial_fractions((
                -ss * p.c,
                ss * p.d,
                0.0,
            ), denominator)
            return (even=even, odd=odd)
    end
    function pole(h, r)
        z = z_of_r(r)
        data = pole_data(h)
        return (_four_complex_even_primitive(data.even, z, p.m) +
            _four_complex_odd_primitive(data.odd, z, p.m)) / p.omega
    end
    function inverse_i0(target)
        z = Elliptic.Jacobi.sn(p.omega * float(target), p.m)
        return radius_from_z(z)
    end
    return (
        kind=:yang_wang_real_jacobi_four_complex,
        roots=(u=p.u, w=p.w, v=p.v, s=p.s),
        parameters=(lambda1=p.lambda1, lambda2=p.lambda2,
            modulus=p.m, alpha=p.alpha, n=p.n, omega=p.omega),
        basis=basis, pole=pole, inverse_i0=inverse_i0,
        infinity_i0=Elliptic.F(asin(p.zinf), p.m) / p.omega,
    )
end

function _c5_polar_sector(a, energy, lz, q, requested)
    candidates = kerr_polar_sector_candidates(a, energy, lz, q)
    if requested === :vortical && candidates == (:constant_latitude,)
        return :constant_latitude
    elseif requested !== nothing
        requested in candidates || error(
            "Requested C5 polar sector $(requested) is not among $(collect(candidates)).")
        return requested
    end
    length(candidates) == 1 || error(
        "C5 polar initial data are ambiguous among $(collect(candidates)).")
    return only(candidates)
end

"""
    kerr_geo_capture_four_complex(a, E, Lz, Q; polar_sector=nothing, polar_hemisphere=:north,
                                  polar_phase=0.0, phi0=0.0, reference_radius=nothing)

The C5 capture (E > 1, four complex radial roots): from infinity into the future horizon at
λ = 0, where τ = v = ψ = 0; t and φ are zero at `reference_radius` (default `λ = λ_∞/2`). On the
spin axis (`polar_sector=:axis_constant`) it is the axis-infall member with azimuth `phi0`.
"""
function kerr_geo_capture_four_complex(
        a::Real, energy::Real, lz::Real, q::Real;
        polar_sector=nothing,
        polar_hemisphere::Symbol=:north,
        polar_phase::Real=0.0,
        phi0::Real=0.0,
        reference_radius=nothing)
    abs(a) < 1 || error(
        "At |a| = 1 the four-complex-root capture is built by `kerr_geo_extremal`.")
    abs(a) > 1.0e-14 || error(
        "Four complex radial roots are absent in the Schwarzschild limit.")
    sector = _c5_polar_sector(a, energy, lz, q, polar_sector)
    sector === :axis_constant && return kerr_geo_capture_axis_infall(a, energy;
        axis=polar_hemisphere, phi0=phi0, reference_radius=reference_radius)
    sector in (:vortical, :constant_latitude, :axis_crossing) || error(
        "C5 supports vortical, constant-latitude, axis-crossing, or exact-axis polar motion.")
    classification = kerr_geo_classify(a, energy, lz, q; polar_sector=sector)
    component = _class_component(classification, :capture, :C5)
    parameters = _four_complex_parameters(energy, classification.Status.root_structure)
    model = _four_complex_model(parameters)
    polar = _polar_solution(a, energy, lz, q, sector, float(polar_phase);
        hemisphere=polar_hemisphere)
    return _direct_component(a, energy, lz, q, component, model, polar, reference_radius;
        roots=(model.roots..., jacobi=model.parameters),
        status=(radial_parameter_residual=parameters.relation_residual,))
end
