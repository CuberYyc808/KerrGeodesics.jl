# Motion along the spin axis (Lz = 0, Q = a²(1 − E²)): the axis-infall members of the Plunge
# (E < 1) and Capture (E ≥ 1) classes, with their shared radial models and endpoint series.

function _axis_legendre_f(phi, m)
    if m <= 1
        return Elliptic.F(phi, m)
    end
    transformed = asin(clamp(sqrt(m) * sin(phi), -1.0, 1.0))
    return Elliptic.F(transformed, inv(m)) / sqrt(m)
end

function _axis_sin_from_f(u, m)
    if m <= 1
        return Elliptic.Jacobi.sn(u, m)
    end
    return Elliptic.Jacobi.sn(sqrt(m) * u, inv(m)) / sqrt(m)
end

function _axis_legendre_basis(phi, m, n; second_order=false)
    if m <= 1
        f = Elliptic.F(phi, m)
        pin = _pi_real(n, phi, m)
        j2 = second_order ? _j2_legendre(n, m, phi) : NaN
        return (F=f, Pi=pin, J2=j2)
    end
    rootm = sqrt(m)
    transformed = asin(clamp(rootm * sin(phi), -1.0, 1.0))
    mt = inv(m)
    nt = n / m
    f = _axis_legendre_f(phi, m)
    pin = _pi_real(nt, transformed, mt) / rootm
    j2 = second_order ?
        _j2_legendre(nt, mt, transformed) / rootm : NaN
    return (F=f, Pi=pin, J2=j2)
end

function _odd_h1(n, m, s)
    u = sqrt(max(1 - m * s^2, 0.0))
    a0 = m - n
    scale = max(1.0, abs(m), abs(n))
    if abs(a0) <= 32 * eps(Float64) * scale
        u > 0 || return Inf
        return 1 / (n * u)
    end
    return -_quadratic_denominator_primitive(a0, n, u)
end

function _odd_h2(n, m, s)
    u = sqrt(max(1 - m * s^2, 0.0))
    a0 = m - n
    scale = max(1.0, abs(m), abs(n))
    if abs(a0) <= 32 * eps(Float64) * scale
        u > 0 || return Inf
        return m / (3 * n^2 * u^3)
    end
    quadratic = _quadratic_denominator_primitive(a0, n, u)
    return -m * (u / (2 * a0 * (a0 + n * u^2)) +
        quadratic / (2 * a0))
end

function _axis_kerr_radial_model(a, energy)
    spin = abs(float(a))
    spin > 0 || error("The Kerr-axis Legendre model requires nonzero spin.")
    k = _e2m1(energy)
    denominator = 1 + spin * k
    denominator > 0 || error("Axis Legendre scale is not real.")
    omega = sqrt(spin * denominator)
    m = 2 / denominator
    chi(r) = atan(r / spin) - pi / 4
    function basis(r)
        amplitude = chi(r)
        s = sin(amplitude)
        legendre2 = _axis_legendre_basis(amplitude, m, 2.0;
            second_order=true)
        h1 = _odd_h1(2.0, m, s)
        h2 = _odd_h2(2.0, m, s)
        return (
            I0=legendre2.F / omega,
            I1=spin * (legendre2.Pi + 2 * h1) / omega,
            I2=spin^2 * (-legendre2.F + 2 * legendre2.J2 +
                4 * h2) / omega,
        )
    end
    function pole(h, r)
        amplitude = chi(r)
        s = sin(amplitude)
        n = 2 * (spin^2 + h^2) / (spin - h)^2
        legendre = _axis_legendre_basis(amplitude, m, n)
        h1 = _odd_h1(n, m, s)
        c0 = -2 * h / n
        c1 = spin - h + 2 * h / n
        return (c0 * legendre.F + c1 * legendre.Pi - 2 * spin * h1) /
            ((spin - h)^2 * omega)
    end
    lambda_primitive(r) = _axis_legendre_f(chi(r), m) / omega
    function radius_from_primitive(value)
        s = clamp(_axis_sin_from_f(omega * value, m), -1.0, 1.0)
        return spin * tan(asin(s) + pi / 4)
    end
    return (
        kind=:axis_legendre_reduction,
        spin=spin,
        k=k,
        m=m,
        omega=omega,
        chi=chi,
        basis=basis,
        pole=pole,
        lambda_primitive=lambda_primitive,
        radius_from_primitive=radius_from_primitive,
    )
end

function _axis_schwarzschild_j1(k, w)
    if k > 0
        root = sqrt(k)
        return log(abs((w - root) / (w + root))) / (2 * root)
    elseif k < 0
        root = sqrt(-k)
        return atan(w / root) / root
    end
    return -inv(w)
end

function _axis_schwarzschild_j2(k, w)
    abs(k) <= 1.0e-14 && return -inv(3 * w^3)
    return -w / (2 * k * (w^2 - k)) -
        _axis_schwarzschild_j1(k, w) / (2 * k)
end

function _axis_schwarzschild_time_primitive(energy, w)
    k = _e2m1(energy)
    return 2 * log(abs((energy + w) / (energy - w))) +
        4 * energy * (_axis_schwarzschild_j1(k, w) +
            _axis_schwarzschild_j2(k, w))
end

function _axis_schwarzschild_v(energy, w)
    k = _e2m1(energy)
    denominator = w^2 - k
    primitive = 4 * log(energy + w) - 2 * log(denominator) +
        4 * energy * (_axis_schwarzschild_j1(k, w) +
            _axis_schwarzschild_j2(k, w)) + 2 / denominator
    horizon = 4 * log(2 * energy) +
        4 * energy * (_axis_schwarzschild_j1(k, energy) +
            _axis_schwarzschild_j2(k, energy)) + 2
    return primitive - horizon
end

# the largest offset from the horizon (≤ 1e-3) where the series' last term is below 1e-15
function _axis_match_offset(coefficients, order)
    match_offset = 1.0e-3
    last_term = Inf
    while match_offset > 1.0e-12
        last_term = abs(coefficients[order] * match_offset^(order + 1) / (order + 1))
        last_term <= 1.0e-15 && break
        match_offset /= 2
    end
    return match_offset, last_term
end

function _axis_endpoint_series(a, energy, q, time_increment)
    horizons = kerr_horizons(a)
    marker = (rplus=horizons.rplus,)
    order = 10
    coefficients = _c3_regular_endpoint_coefficients(
        a, energy, 0.0, q, marker, :v, order)
    match_offset, last_term = _axis_match_offset(coefficients, order)
    last_term <= 1.0e-15 || error(
        "Axis endpoint series did not reach the machine-precision match target.")
    rmatch = horizons.rplus + match_offset
    local_value = _c3_regular_endpoint_series_integral(
        coefficients, match_offset, order)
    rstar(r) = kerr_rstar(a, r)
    function endpoint(r)
        y = r - horizons.rplus
        y >= -1.0e-13 || throw(DomainError(
            r, "Axis endpoint radius lies inside the future horizon."))
        y <= 1.0e-15 && return 0.0
        if y <= match_offset
            return _c3_regular_endpoint_series_integral(
                coefficients, y, order)
        end
        return local_value + time_increment(rmatch, r) -
            (rstar(r) - rstar(rmatch))
    end
    return (endpoint=endpoint, order=order, match_offset=match_offset,
        last_term_abs=last_term, coefficients=coefficients)
end

# Radial pieces of an axis infall: Schwarzschild elementary forms or the Kerr Legendre model.
# Each branch is its own function so every captured variable is assigned once.
_axis_radial_parts(spin, evalue, qaxis, rplus) = abs(spin) <= 1.0e-14 ?
    _axis_schwarzschild_parts(evalue) :
    _axis_kerr_parts(_axis_kerr_radial_model(spin, evalue), spin, evalue, qaxis, rplus)

function _axis_schwarzschild_parts(evalue)
    k = evalue^2 - 1
    w = lambda -> evalue + float(lambda)
    return (model=nothing,
        lambda_start=evalue >= 1 ? sqrt(max(k, 0.0)) - evalue : -evalue,
        radius_from_lambda=lambda -> 2 / (w(lambda)^2 - k),
        start_radius=evalue < 1 ? 2 / (1 - evalue^2) : Inf,
        time_increment=(left, right) ->
            _axis_schwarzschild_time_primitive(evalue, sqrt(k + 2 / left)) -
            _axis_schwarzschild_time_primitive(evalue, sqrt(k + 2 / right)),
        proper_increment=(left, right) -> 4 * (
            _axis_schwarzschild_j2(k, sqrt(k + 2 / left)) -
            _axis_schwarzschild_j2(k, sqrt(k + 2 / right))),
        regular_v_endpoint=r -> -_axis_schwarzschild_v(evalue, sqrt(k + 2 / r)),
        formula_kind=:schwarzschild_axis_elementary,
        endpoint_series=nothing)
end

function _axis_kerr_parts(model, spin, evalue, qaxis, rplus)
    lambda_horizon_primitive = model.lambda_primitive(rplus)
    deficit = 1 - evalue^2
    start_radius = evalue < 1 ? (1 + sqrt(1 - deficit^2 * spin^2)) / deficit : Inf
    lambda_start = lambda_horizon_primitive - model.lambda_primitive(start_radius)
    radius_from_lambda = function(lambda)
        lam = float(lambda)
        abs(lam) <= 2.0e-14 && return rplus
        evalue < 1 && abs(lam - lambda_start) <= 2.0e-13 &&
            return start_radius
        return model.radius_from_primitive(lambda_horizon_primitive - lam)
    end
    increments = evalue == 1 ? _axis_parabolic_increments(spin, qaxis) :
        _axis_legendre_increments(model, spin, evalue)
    endpoint_series = _axis_endpoint_series(spin, evalue, qaxis, increments.time)
    return (; model, lambda_start, radius_from_lambda, start_radius,
        time_increment=increments.time, proper_increment=increments.proper,
        regular_v_endpoint=endpoint_series.endpoint, formula_kind=increments.kind,
        endpoint_series)
end

function _axis_parabolic_increments(spin, qaxis)
    cparams = _c1_one_real_parameters(spin, 0.0, qaxis)
    cparams === nothing && error(
        "The E = 1 axis radial cubic must have one real root below r+ and a complex pair.")
    return (time=(left, right) -> _c1_radial_time_increment(spin, 0.0, cparams, left, right),
        proper=(left, right) -> 0.5 * _c1_k_legendre_delta(cparams, 2, left, right) +
            2 * spin^2 * _c1_k_legendre_delta(cparams, 0, left, right),
        kind=:axis_parabolic_cubic_legendre)
end

function _axis_legendre_increments(model, spin, evalue)
    residues = _radial_residues(spin, evalue, 0.0)
    function time_primitive(r)
        basis = model.basis(r)
        return evalue * basis.I2 + 2 * evalue * basis.I1 +
            (evalue * spin^2 + 4 * evalue) * basis.I0 +
            residues.c_t_plus * model.pole(residues.rplus, r) +
            residues.c_t_minus * model.pole(residues.rminus, r)
    end
    function proper(left, right)
        left_basis = model.basis(left)
        right_basis = model.basis(right)
        return (right_basis.I2 + spin^2 * right_basis.I0) -
            (left_basis.I2 + spin^2 * left_basis.I0)
    end
    return (time=(left, right) -> time_primitive(right) - time_primitive(left),
        proper=proper, kind=evalue < 1 ? :axis_elliptic_legendre : :axis_hyperbolic_legendre)
end

function _axis_infall_member(a::Real, energy::Real;
        axis,
        broad_class,
        phi0::Real=0.0,
        reference_radius=nothing)
    axis in (:north, :south) || error("Specify axis=:north or axis=:south.")
    0 < energy || error("Future-directed axis infall requires E>0.")
    abs(a) < 1 || error("This constructor requires |a| < 1; axis infall at |a| = 1 is " *
        "built by kerr_geo_extremal(a, E, 0, a^2(1 - E^2); axis=axis).")
    broad_class === :plunge && energy < 1 ||
        broad_class === :capture && energy >= 1 || error(
            "Plunge axis infall requires E<1; capture axis infall requires E>=1.")

    spin = float(a)
    evalue = float(energy)
    qaxis = kerr_axis_carter_q(spin, evalue)
    horizons = kerr_horizons(spin)
    rplus = horizons.rplus
    classification = kerr_geo_classify(
        spin, evalue, 0.0, qaxis; polar_sector=:axis_constant)
    component = only(c for c in classification.Components if c.BroadClass === broad_class)

    (; model, lambda_start, radius_from_lambda, start_radius, time_increment,
        proper_increment, regular_v_endpoint, formula_kind, endpoint_series) =
        _axis_radial_parts(spin, evalue, qaxis, rplus)

    function check_regular(lambda)
        lam = float(lambda)
        lambda_start < lam <= 0 ||
            (evalue < 1 && isapprox(lam, lambda_start; atol=2.0e-13)) ||
            throw(DomainError(lambda,
                "Mino time must lie between λ = $(lambda_start) and the future horizon λ = 0."))
        return clamp(lam, lambda_start, 0.0)
    end
    function check_bl(lambda)
        lam = check_regular(lambda)
        valid_start = evalue < 1 ? lambda_start <= lam : lambda_start < lam
        valid_start && lam < 0 || throw(DomainError(
            lambda, "BL t and r* exclude the future-horizon endpoint λ = 0."))
        return lam
    end
    r(lambda) = radius_from_lambda(check_regular(lambda))
    theta_value = axis === :north ? 0.0 : pi
    theta(lambda) = (check_regular(lambda); theta_value)
    z(lambda) = (check_regular(lambda); axis === :north ? 1.0 : -1.0)
    phi(lambda) = (check_regular(lambda); float(phi0))
    lambda_reference = if reference_radius === nothing
        evalue < 1 ? lambda_start : 0.5 * lambda_start
    else
        rref_input = float(reference_radius)
        rref_input > rplus || error("Axis BL reference radius must be exterior.")
        if model === nothing
            sqrt(evalue^2 - 1 + 2 / rref_input) - evalue
        else
            model.lambda_primitive(rplus) - model.lambda_primitive(rref_input)
        end
    end
    rref = isfinite(start_radius) && lambda_reference == lambda_start ?
        start_radius : r(lambda_reference)
    t(lambda) = time_increment(r(check_bl(lambda)), rref)
    v(lambda) = -regular_v_endpoint(r(check_regular(lambda)))
    tau(lambda) = -proper_increment(rplus, r(check_regular(lambda)))
    rstar(lambda) = kerr_rstar(spin, r(check_bl(lambda)))

    radial_potential(rvalue) = kerr_axis_radial_potential(spin, evalue, rvalue)
    polar_potential(zvalue) = kerr_polar_z_potential(spin, evalue, 0.0, qaxis, zvalue)
    ur(lambda) = -sqrt(max(radial_potential(r(lambda)), 0.0))
    vanishing(lambda) = (check_regular(lambda); 0.0)
    # dt/dλ = E Σ²/Δ on the axis (Σ = r² + a²), and the proper-time velocity for the norm
    ut(lambda) = (rv = r(check_bl(lambda)); evalue * (rv^2 + spin^2)^2 / kerr_delta(spin, rv))
    function normalization_residual(lambda)
        rv = r(check_bl(lambda))
        sigma = rv^2 + spin^2
        delta = kerr_delta(spin, rv)
        return -delta / sigma * (ut(lambda) / sigma)^2 + sigma / delta * (ur(lambda) / sigma)^2 + 1
    end
    start_role = evalue < 1 ? :finite_turning_point : :past_infinity
    return _member(broad_class, component.CaseId; component=component,
        constants=(a=spin, E=evalue, Lz=0.0, Q=qaxis),
        roots=(radial=Tuple(item.radius for item in classification.Status.root_structure.real_roots),
            polar=(sector=:axis_constant, axis=axis)),
        reference=(lambda0_event=:future_horizon, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=lambda_reference, t_phi_zero_radius=rref,
            tau_zero_event=:future_horizon, lambda_regular=0.0, phi0=float(phi0)),
        domain=(mino=(lambda_start, 0.0), endpoint_closed=(evalue < 1, true),
            endpoint_roles=(start_role, :future_horizon), horizon_lambda=0.0),
        # on the axis φ is a gauge (phi0) and so is ψ: ψ = φ, not φ + φ_H
        trajectory=(t=t, r=r, theta=theta, z=z, phi=phi, tau=tau, rstar=rstar, v=v, psi=phi),
        velocity=(ut=ut, ur=ur, uz=vanishing, utheta=vanishing, uphi=vanishing,
            dtau_dlambda=lambda -> r(check_regular(lambda))^2 + spin^2),
        potentials=(radial=radial_potential, polar_z=polar_potential),
        residuals=(radial=lambda -> ur(lambda)^2 - radial_potential(r(lambda)),
            polar_z=lambda -> polar_potential(z(lambda)),
            normalization=normalization_residual),
        status=(supported=true, name=evalue < 1 ? :finite_axis_infall : :infinity_axis_infall,
            formula_family=component.FormulaFamily, formula_kind=formula_kind, q_axis=qaxis,
            endpoint_series_order=endpoint_series === nothing ? 0 : endpoint_series.order,
            endpoint_series_match_offset=endpoint_series === nothing ?
                0.0 : endpoint_series.match_offset,
            endpoint_series_last_term_abs=endpoint_series === nothing ?
                0.0 : endpoint_series.last_term_abs,
            polar=(sector=:axis_constant, axis=axis)))
end

"""
    kerr_geo_plunge_axis_infall(a, E; axis, phi0=0.0, reference_radius=nothing)

Motion along the spin axis (`axis = :north` or `:south`; Lz = 0, Q = a²(1 − E²)) from the
turning point into the future horizon, for 0 < E < 1 and |a| < 1. λ = 0 on the future horizon,
where τ = v = 0; t vanishes at `reference_radius` (default: the turning point); φ = ψ = `phi0`.
"""
function kerr_geo_plunge_axis_infall(a::Real, energy::Real;
        axis,
        phi0::Real=0.0,
        reference_radius=nothing)
    return _axis_infall_member(
        a, energy; axis=axis,
        broad_class=:plunge, phi0=phi0,
        reference_radius=reference_radius)
end

"""
    kerr_geo_capture_axis_infall(a, E; axis, phi0=0.0, reference_radius=nothing)

Motion along the spin axis (`axis = :north` or `:south`; Lz = 0, Q = a²(1 − E²)) from infinity
into the future horizon, for E ≥ 1 and |a| < 1. λ = 0 on the future horizon, where τ = v = 0;
t vanishes at `reference_radius` (default: the radius halfway in Mino time between infinity
and the horizon); φ = ψ = `phi0`.
"""
kerr_geo_capture_axis_infall(a::Real, energy::Real; axis, phi0::Real=0.0,
        reference_radius=nothing) = _axis_infall_member(a, energy; axis=axis,
    broad_class=:capture, phi0=phi0, reference_radius=reference_radius)
