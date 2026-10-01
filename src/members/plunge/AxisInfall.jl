# Motion along the spin axis (Lz = 0, Q = a²(1 − E²)): the axis-infall members of the Plunge
# (E < 1) and Capture (E ≥ 1) classes. r(λ) comes from the elementary (a = 0) or Legendre (a ≠ 0)
# radial model; t, τ, v from the radial engine, as for every other member. On the axis φ is a
# gauge (phi0), and so is ψ.

# Legendre forms of the axis model, whose parameter m = 2/(1 + a(E² − 1)) exceeds 1 for
# a(E² − 1) < 1: there the reciprocal-modulus transformation F(φ|m) = F(ψ|1/m)/√m,
# sin ψ = √m sin φ, applies, with the complement 1 − 1/m = −m1/m
function _axis_legendre_f(phi, m, m1)
    m <= 1 && return _ellip_f(phi, m1)
    return _ellip_f(asin(clamp(sqrt(m) * sin(phi), -1.0, 1.0)), -m1 / m) / sqrt(m)
end

# the Landen record of the parameter actually used, m or 1/m
_axis_landen(m, m1) = m <= 1 ? _landen(m, m1) : _landen(inv(m), -m1 / m)

_axis_sin_from_f(u, L, m) = m <= 1 ? _ellipj_reduced(u, L)[1] :
    _ellipj_reduced(sqrt(m) * u, L)[1] / sqrt(m)

function _axis_legendre_basis(phi, m, m1, n; second_order=false)
    if m <= 1
        f = _ellip_f(phi, m1)
        pin = _ellip_pi(phi, m1, n, 1 - n)
        j2 = second_order ? _ellip_pi2(phi, m1, n, 1 - n) : NaN
        return (F=f, Pi=pin, J2=j2)
    end
    rootm = sqrt(m)
    transformed = asin(clamp(rootm * sin(phi), -1.0, 1.0))
    mt, mt1 = inv(m), -m1 / m
    nt = n / m
    f = _axis_legendre_f(phi, m, m1)
    pin = _ellip_pi(transformed, mt1, nt, 1 - nt) / rootm
    j2 = second_order ? _ellip_pi2(transformed, mt1, nt, 1 - nt) / rootm : NaN
    return (F=f, Pi=pin, J2=j2)
end

# The odd-power primitives ∫ ds/((1 − n s²)^k √(1 − m s²)) in u = √(1 − m s²), whose quadratic
# denominator a0 + n u² has a0 = m − n: within rounding of a0 = 0 the primitive is the
# elementary limit (the atan form's constant π/(2√(a0 n)) diverges there).
function _odd_argument(n, m, s)
    u = sqrt(max(1 - m * s^2, 0.0))
    a0 = m - n
    return u, a0, abs(a0) <= 32 * eps(Float64) * max(1.0, abs(m), abs(n))
end

function _odd_h1(n, m, s)
    u, a0, limit = _odd_argument(n, m, s)
    if limit
        u > 0 || return Inf
        return 1 / (n * u)
    end
    return -_quadratic_denominator_primitive(a0, n, u)
end

function _odd_h2(n, m, s)
    u, a0, limit = _odd_argument(n, m, s)
    if limit
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
    m1 = (spin * k - 1) / denominator
    L = _axis_landen(m, m1)
    chi(r) = atan(r / spin) - pi / 4
    function basis(r)
        amplitude = chi(r)
        s = sin(amplitude)
        legendre2 = _axis_legendre_basis(amplitude, m, m1, 2.0;
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
        legendre = _axis_legendre_basis(amplitude, m, m1, n)
        h1 = _odd_h1(n, m, s)
        c0 = -2 * h / n
        c1 = spin - h + 2 * h / n
        return (c0 * legendre.F + c1 * legendre.Pi - 2 * spin * h1) /
            ((spin - h)^2 * omega)
    end
    lambda_primitive(r) = _axis_legendre_f(chi(r), m, m1) / omega
    function radius_from_primitive(value)
        s = clamp(_axis_sin_from_f(omega * value, L, m), -1.0, 1.0)
        return spin * tan(asin(s) + pi / 4)
    end
    return (
        kind=:axis_legendre_reduction,
        spin=spin,
        k=k,
        m=m,
        m1=m1,
        omega=omega,
        chi=chi,
        basis=basis,
        pole=pole,
        lambda_primitive=lambda_primitive,
        radius_from_primitive=radius_from_primitive,
    )
end







# Radial pieces of an axis infall: the Mino time from the start (turning point or infinity) to
# the horizon, r at the Mino time δ before the horizon, its inverse, and the start radius.
_axis_radial_parts(spin, evalue, rplus) = kerr_metric_limit(spin) === :schwarzschild ?
    _axis_schwarzschild_parts(evalue) :
    _axis_kerr_parts(_axis_kerr_radial_model(spin, evalue), evalue, rplus)

function _axis_schwarzschild_parts(evalue)
    k = _e2m1(evalue)
    w = lambda -> evalue + float(lambda)
    # δ is the Mino time before the horizon: r = 2/((E − δ)² − k), δ = E − √(k + 2/r)
    return (lambda_start=evalue >= 1 ? sqrt(max(k, 0.0)) - evalue : -evalue,
        radius=δ -> 2 / (w(-δ)^2 - k),
        mino=r -> evalue - sqrt(max(k + 2 / r, 0.0)),
        start_radius=evalue < 1 ? -2 / k : Inf,
        formula_kind=:schwarzschild_axis_elementary)
end

function _axis_kerr_parts(model, evalue, rplus)
    lambda_horizon_primitive = model.lambda_primitive(rplus)
    deficit = -_e2m1(evalue)
    start_radius = evalue < 1 ? (1 + sqrt(1 - deficit^2 * model.spin^2)) / deficit : Inf
    lambda_start = lambda_horizon_primitive - model.lambda_primitive(start_radius)
    # δ is the Mino time before the horizon (δ = 0 there)
    radius = function(δ)
        δ <= 0 && return rplus
        evalue < 1 && δ >= -lambda_start && return start_radius
        return clamp(model.radius_from_primitive(lambda_horizon_primitive + δ), rplus,
            start_radius)
    end
    mino(r) = model.lambda_primitive(r) - lambda_horizon_primitive
    return (; lambda_start, radius, mino, start_radius, formula_kind=:axis_legendre)
end

# The classified component of class `broad_class` of axis constants (Lz = 0, Q = a²(1 − E²)):
# from a classification (the family passes its own), or classified here for the constructors.
function _axis_component(classification, broad_class)
    components = [c for c in classification.Components if c.BroadClass === broad_class]
    isempty(components) && error("These constants admit no $(kerr_geo_class(broad_class).name) " *
        "member on the spin axis; their cases are $(classification.CaseIds).")
    return only(components)
end

function _axis_component(a, energy, broad_class)
    abs(a) < 1 || error("This constructor requires |a| < 1; axis infall at |a| = 1 is " *
        "built by kerr_geo_extremal(a, E, 0, a^2(1 - E^2); axis=axis).")
    broad_class === :plunge && energy < 1 ||
        broad_class === :capture && energy >= 1 || error(
            "Plunge axis infall requires E<1; capture axis infall requires E>=1.")
    spin, evalue = float(a), float(energy)
    classification = kerr_geo_classify(spin, evalue, 0.0, kerr_axis_carter_q(spin, evalue);
        polar_sector=:axis_constant)
    return _axis_component(classification, broad_class)
end

function _axis_infall_member(a::Real, energy::Real, component;
        axis,
        phi0::Real=0.0,
        reference_radius=nothing)
    axis in (:north, :south) || error("Specify axis=:north or axis=:south.")
    broad_class = component.BroadClass

    spin = float(a)
    evalue = float(energy)
    qaxis = kerr_axis_carter_q(spin, evalue)
    rplus = kerr_horizons(spin).rplus
    (; lambda_start, radius, mino, start_radius, formula_kind) =
        _axis_radial_parts(spin, evalue, rplus)

    function check_regular(lambda)
        lam = float(lambda)
        lambda_start < lam <= 0 ||
            (evalue < 1 && abs(lam - lambda_start) <= MINO_ENDPOINT_TOL) ||
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
    r_of(lambda) = radius(-lambda)                   # δ = −λ is the Mino time before the horizon
    r(lambda) = r_of(check_regular(lambda))
    z0 = axis === :north ? 1.0 : -1.0
    theta(lambda) = (check_regular(lambda); acos(z0))
    z(lambda) = (check_regular(lambda); z0)
    phi(lambda) = (check_regular(lambda); float(phi0))
    lambda_reference = if reference_radius === nothing
        evalue < 1 ? lambda_start : 0.5 * lambda_start
    else
        rref_input = float(reference_radius)
        rref_input > rplus || error("Axis BL reference radius must be exterior.")
        -mino(rref_input)
    end
    rref = isfinite(start_radius) && lambda_reference == lambda_start ?
        start_radius : r(lambda_reference)
    # the polar motion is the constant z = ±1: its t and φ rates vanish (Lz = 0), dτ/dλ = a²
    polar = (formula=lambda -> (z=z0, uz=0.0, sin2=0.0, theta=acos(z0), phi=0.0, t=0.0,
            tau=spin^2 * float(lambda)),
        metadata=(sector=:axis_constant, axis=axis))
    coords = _engine_coordinates(spin, evalue, 0.0, qaxis, r_of, _polar_primitive(polar);
        potential=_radial_potential_from_roots(spin, evalue, 0.0, qaxis, component.Metadata.structure),
        domain=(lambda_start, 0.0), ends=(evalue < 1 ? :turning : :infinity, :horizon), σ=-1.0,
        λ_bl=lambda_reference, λ_regular=0.0, σ_regular=-1.0)
    t(lambda) = _coords_t(coords, check_bl(lambda))
    v(lambda) = _coords_v(coords, check_regular(lambda))
    tau(lambda) = _coords_tau(coords, check_regular(lambda))
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
    return _member(broad_class, component.CaseId, kerr_geo_tier(component.CaseId), component,
        (a=spin, E=evalue, Lz=0.0, Q=qaxis),
        (radial=Tuple(item.radius for item in component.Metadata.roots),
            polar=polar.metadata),
        (lambda0_event=:future_horizon, t_phi_zero_event=:reference_radius,
            t_phi_zero_lambda=lambda_reference, t_phi_zero_radius=rref,
            tau_zero_event=:future_horizon, lambda_regular=0.0, phi0=float(phi0)),
        (mino=(lambda_start, 0.0), endpoint_closed=(evalue < 1, true),
            endpoint_roles=(start_role, :future_horizon), horizon_lambda=0.0),
        (t=t, r=r, theta=theta, z=z, phi=phi, tau=tau, rstar=rstar, v=v, psi=phi),
        (ut=ut, ur=ur, uz=vanishing, utheta=vanishing, uphi=vanishing,
            dtau_dlambda=lambda -> r(check_regular(lambda))^2 + spin^2),
        (radial=radial_potential, polar_z=polar_potential),
        (radial=lambda -> ur(lambda)^2 - radial_potential(r(lambda)),
            polar_z=lambda -> polar_potential(z(lambda)),
            normalization=normalization_residual),
        (supported=true, name=evalue < 1 ? :finite_axis_infall : :infinity_axis_infall,
            formula_family=component.FormulaFamily, formula_kind=formula_kind, q_axis=qaxis,
            polar=polar.metadata),
        SpectralStatus(() -> (_coords_spectral(coords), _polar_spectral(polar))))
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
    return _axis_infall_member(a, energy, _axis_component(a, energy, :plunge); axis=axis,
        phi0=phi0, reference_radius=reference_radius)
end

"""
    kerr_geo_capture_axis_infall(a, E; axis, phi0=0.0, reference_radius=nothing)

Motion along the spin axis (`axis = :north` or `:south`; Lz = 0, Q = a²(1 − E²)) from infinity
into the future horizon, for E ≥ 1 and |a| < 1. λ = 0 on the future horizon, where τ = v = 0;
t vanishes at `reference_radius` (default: the radius halfway in Mino time between infinity
and the horizon); φ = ψ = `phi0`.
"""
kerr_geo_capture_axis_infall(a::Real, energy::Real; axis, phi0::Real=0.0,
        reference_radius=nothing) = _axis_infall_member(a, energy,
    _axis_component(a, energy, :capture); axis=axis, phi0=phi0, reference_radius=reference_radius)
