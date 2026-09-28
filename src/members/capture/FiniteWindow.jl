# Finite-window capture API `kerr_geo_capture` (E >= 1, infinity -> future horizon) and the
# closed-form r(λ), λ(r) of the E = 1 (C1) and E > 1 (C3) captures it uses.

"""
    KerrGeoCapture

E ≥ 1 orbit falling from infinity into the future horizon. `Formula` is
`:hyperbolic_capture` (E > 1) or `:parabolic_capture` (E = 1). Mino time λ is zero on
the horizon; the trajectory gives r, θ, the regular horizon coordinates v and ψ, and
finite-window increments of the Boyer–Lindquist t and φ (which diverge on the horizon).
"""
struct KerrGeoCapture
    Formula::Symbol
    EnergyRegime::Symbol
    Outcome::Symbol
    OrbitalParameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    Parametrization::String
    Roots::NamedTuple
    ReferenceZero::NamedTuple
    Trajectory::Any
    Velocity::Any
    Potentials::NamedTuple
    Residuals::NamedTuple
    Status::NamedTuple
end

function _capture_domain(lambda_infinity)
    return (mino=(lambda_infinity, 0.0), endpoint_roles=(:past_infinity, :future_horizon),
            endpoint_closed=(false, true))
end

function _c3_complex_parameters(a, energy, lz, q; atol=1e-10)
    roots_all = radial_roots_for_constants(a, energy, lz, q)
    real_roots = Float64[]
    complex_roots = ComplexF64[]
    for root in roots_all
        if abs(imag(root)) <= atol
            push!(real_roots, real(root))
        else
            push!(complex_roots, root)
        end
    end
    sort!(real_roots)
    length(real_roots) == 2 || return nothing
    length(complex_roots) == 2 || return nothing
    upper = complex_roots[argmax(imag.(complex_roots))]
    eta = abs(imag(upper))
    eta > 0 || return nothing
    rplus = _rplus(a)
    real_roots[1] < real_roots[2] < rplus || return nothing
    lead = _e2m1(energy)
    lead > 0 || return nothing
    return (
        lead=lead,
        r1=real_roots[1],
        r2=real_roots[2],
        rho=real(upper),
        eta=eta,
        rplus=rplus,
        roots_all=Tuple(roots_all),
        real_roots=Tuple(real_roots),
    )
end

function _c3_shape(c)
    b = sqrt((c.r2 - c.rho)^2 + c.eta^2)
    d = sqrt((c.r1 - c.rho)^2 + c.eta^2)
    n = ((c.r2 - c.rho) * (c.rho - c.r1) - c.eta^2) / (b * d)
    m = (1 - n) / 2
    return (B=b, C=d, n=n, m=m)
end

function _c3_phi_of_r(c, r)
    shape = _c3_shape(c)
    s2 = (r - c.r2) / (r - c.r1)
    y = sqrt(shape.C / shape.B) * sqrt(max(s2, 0.0))
    return 2 * atan(y)
end

function _c3_lambda_of_r(c, r)
    shape = _c3_shape(c)
    return Elliptic.F(_c3_phi_of_r(c, r), shape.m) /
           (sqrt(c.lead) * sqrt(shape.B * shape.C))
end

function _c3_lambda_infinity(c)
    shape = _c3_shape(c)
    phi_infinity = 2 * atan(sqrt(shape.C / shape.B))
    return Elliptic.F(phi_infinity, shape.m) /
           (sqrt(c.lead) * sqrt(shape.B * shape.C))
end

function _capture_positive_radial_interval(r_left, r_right)
    r_left == r_right && return (r_left, r_right, 1)
    r_left < r_right && return (r_left, r_right, 1)
    return (r_right, r_left, -1)
end

function _c3_bl_laurent_coefficients(a, energy, lz, q, c, component, order)
    rp = c.rplus
    rm = _rminus(a)
    d = rp - rm
    pplus = 2 * energy * rp - a * lz
    pplus > 0 || error("C3 endpoint anchor requires positive future-horizon P_+.")
    bconst = (lz - a * energy)^2 + q
    inv_order = order + 1

    p = zeros(Float64, inv_order + 1)
    p[1] = pplus
    p[2] = 2 * energy * rp
    p[3] = energy

    delta = zeros(Float64, inv_order + 1)
    delta[2] = d
    delta[3] = 1.0

    sgeom = zeros(Float64, inv_order + 1)
    sgeom[1] = rp^2 + bconst
    sgeom[2] = 2 * rp
    sgeom[3] = 1.0

    rseries = _series_mul(p, p, inv_order) .-
              _series_mul(delta, sgeom, inv_order)
    invsqrt = _series_inv(_series_sqrt_positive(rseries, inv_order), inv_order)

    regular_denom = zeros(Float64, inv_order + 1)
    regular_denom[1] = d
    regular_denom[2] = 1.0

    dot = zeros(Float64, order + 2) # powers -1:order, offset by +2.
    if component === :psi
        p_over_regular = _series_div(p, regular_denom, inv_order)
        dot[1] = a * p_over_regular[1]
        for k in 0:order
            dot[k + 2] = a * p_over_regular[k + 2] -
                         (k == 0 ? a * energy : 0.0)
        end
    elseif component === :v
        r2pa2 = zeros(Float64, inv_order + 1)
        r2pa2[1] = 2 * rp
        r2pa2[2] = 2 * rp
        r2pa2[3] = 1.0
        numerator = _series_mul(r2pa2, p, inv_order)
        t_over_regular = _series_div(numerator, regular_denom, inv_order)
        dot[1] = t_over_regular[1]
        for k in 0:order
            dot[k + 2] = t_over_regular[k + 2]
        end
    else
        error("Unknown C3 endpoint component $(component).")
    end

    coeffs = Dict{Int,Float64}()
    for k in -1:order
        s = 0.0
        for j in -1:k
            inv_idx = k - j
            0 <= inv_idx <= inv_order || continue
            s += dot[j + 2] * invsqrt[inv_idx + 1]
        end
        coeffs[k] = s
    end
    return coeffs
end

function _c3_subtraction_laurent_coefficients(a, c, component, order)
    rp = c.rplus
    rm = _rminus(a)
    d = rp - rm
    denom = zeros(Float64, order + 2)
    denom[1] = d
    denom[2] = 1.0
    coeffs = Dict{Int,Float64}()
    if component === :psi
        q = _series_inv(denom, order + 1)
        coeffs[-1] = a * q[1]
        for k in 0:order
            coeffs[k] = a * q[k + 2]
        end
        return coeffs
    elseif component === :v
        numerator = zeros(Float64, order + 2)
        numerator[1] = 2 * rp
        numerator[2] = 2 * rp
        numerator[3] = 1.0
        q = _series_div(numerator, denom, order + 1)
        coeffs[-1] = q[1]
        for k in 0:order
            coeffs[k] = q[k + 2]
        end
        return coeffs
    end
    error("Unknown C3 endpoint component $(component).")
end

function _c3_regular_endpoint_coefficients(a, energy, lz, q, c, component, order)
    bl = _c3_bl_laurent_coefficients(a, energy, lz, q, c, component, order)
    sub = _c3_subtraction_laurent_coefficients(a, c, component, order)
    abs(bl[-1] - sub[-1]) <= 1.0e-10 ||
        error("C3 regular endpoint pole cancellation failed for $(component).")
    coeffs = Dict{Int,Float64}()
    for k in 0:order
        coeffs[k] = bl[k] - sub[k]
    end
    return coeffs
end

function _c3_regular_endpoint_series_integral(coeffs, y, order)
    value = 0.0
    for k in 0:order
        value += coeffs[k] * y^(k + 1) / (k + 1)
    end
    return value
end

function _c3_radius_from_horizon_lambda(c, lambda; max_iter=90)
    lambda_horizon = _c3_lambda_of_r(c, c.rplus)
    lambda_infinity = lambda_horizon - _c3_lambda_infinity(c)
    lambda <= 2e-13 || error("C3 horizon-zero lambda must be nonpositive outside the future horizon.")
    lambda >= lambda_infinity - 2e-13 ||
        error("C3 lambda is beyond the infinity endpoint for this finite branch.")
    abs(lambda) <= 2e-13 && return c.rplus
    target = lambda_horizon - lambda
    # closed-form inverse of λ(r) = F(φ(r)|m)/sqrt(lead B C): φ = am(...), then
    # tan^2(φ/2) B/C = (r - r2)/(r - r1)
    shape = _c3_shape(c)
    φ = Elliptic.Jacobi.am(target * sqrt(c.lead * shape.B * shape.C), shape.m)
    ratio = tan(φ / 2)^2 * shape.B / shape.C
    closed = (c.r2 - ratio * c.r1) / (1 - ratio)
    isfinite(closed) && closed >= c.rplus && return closed
    low = c.rplus
    high = c.rplus + 1.0
    while _c3_lambda_of_r(c, high) < target
        high *= 1.5
        high > 1e10 && error("Failed to bracket C3 radius from horizon-zero Mino time.")
    end
    for _ in 1:max_iter
        mid = 0.5 * (low + high)
        if _c3_lambda_of_r(c, mid) < target
            low = mid
        else
            high = mid
        end
    end
    return 0.5 * (low + high)
end

function _c1_one_real_parameters(a, lz, q; atol=1e-10)
    roots_all = radial_roots_for_constants(a, 1.0, lz, q)
    real_roots = Float64[]
    complex_roots = ComplexF64[]
    for root in roots_all
        if abs(imag(root)) <= atol
            push!(real_roots, real(root))
        else
            push!(complex_roots, root)
        end
    end
    sort!(real_roots)
    length(real_roots) == 1 || return nothing
    length(complex_roots) == 2 || return nothing
    upper = complex_roots[argmax(imag.(complex_roots))]
    eta = abs(imag(upper))
    eta > 0 || return nothing
    rplus = _rplus(a)
    real_roots[1] < rplus || return nothing
    pplus = rplus^2 + a^2 - a * lz
    pplus > 0 || return nothing
    return (
        x0=real_roots[1],
        rho=real(upper),
        eta=eta,
        rplus=rplus,
        pplus=pplus,
        roots_all=Tuple(roots_all),
            real_roots=Tuple(real_roots),
    )
end

function _c1_shape(c)
    d = c.x0 - c.rho
    B = sqrt(d^2 + c.eta^2)
    m = (B - d) / (2 * B)
    return (d=d, B=B, m=m)
end

function _c1_psi_of_r(c, r)
    r > c.x0 || error("C1 radius must lie above the real cubic root.")
    shape = _c1_shape(c)
    return 2 * atan(sqrt(r - c.x0) / sqrt(shape.B))
end

function _c1_lambda_of_r(c, r)
    shape = _c1_shape(c)
    return Elliptic.F(_c1_psi_of_r(c, r), shape.m) / sqrt(2 * shape.B)
end

function _c1_lambda_infinity(c)
    shape = _c1_shape(c)
    return Elliptic.F(pi, shape.m) / sqrt(2 * shape.B)
end

function _c1_u_of_r(c, r)
    r > c.x0 || error("C1 radius must lie above the real cubic root.")
    shape = _c1_shape(c)
    return sqrt(r - c.x0) / sqrt(shape.B)
end

function _c1_q4(shape, u)
    return u^4 + (2 - 4 * shape.m) * u^2 + 1
end

function _c1_k_legendre_delta(c, power, r_left, r_right)
    shape = _c1_shape(c)
    u_left = _c1_u_of_r(c, r_left)
    u_right = _c1_u_of_r(c, r_right)
    psi_left = _c1_psi_of_r(c, r_left)
    psi_right = _c1_psi_of_r(c, r_right)
    delta_f = Elliptic.F(psi_right, shape.m) - Elliptic.F(psi_left, shape.m)
    k0 = 0.5 * delta_f
    power == 0 && return k0
    boundary(u) = u * sqrt(_c1_q4(shape, u)) / (1 + u^2)
    delta_e = Elliptic.E(psi_right, shape.m) - Elliptic.E(psi_left, shape.m)
    k1 = k0 + (boundary(u_right) - boundary(u_left)) - delta_e
    power == 1 && return k1
    avec = 2 - 4 * shape.m
    radial_boundary(u) = u * sqrt(_c1_q4(shape, u))
    k2 = ((radial_boundary(u_right) - radial_boundary(u_left)) -
          2 * avec * k1 - k0) / 3
    power == 2 && return k2
    error("C1 K_j is defined for the powers 0, 1 and 2.")
end

function _c1_radius_from_horizon_lambda(c, lambda; max_iter=90)
    lambda_horizon = _c1_lambda_of_r(c, c.rplus)
    lambda_infinity = lambda_horizon - _c1_lambda_infinity(c)
    lambda <= 2e-13 || error("C1 horizon-zero lambda must be nonpositive outside the future horizon.")
    lambda >= lambda_infinity - 2e-13 ||
        error("C1 lambda is beyond the infinity endpoint for this finite branch.")
    abs(lambda) <= 2e-13 && return c.rplus
    target = lambda_horizon - lambda
    # closed-form inverse of λ(r) = F(ψ(r)|m)/sqrt(2B): r = x0 + B tan^2(ψ/2)
    shape = _c1_shape(c)
    ψ = Elliptic.Jacobi.am(target * sqrt(2 * shape.B), shape.m)
    closed = c.x0 + shape.B * tan(ψ / 2)^2
    isfinite(closed) && closed >= c.rplus && return closed
    low = c.rplus
    span = 1.0
    high = c.rplus + span
    while _c1_lambda_of_r(c, high) < target
        span *= 1.5
        high = c.rplus + span
        high > 1e10 && error("Failed to bracket C1 radius from horizon-zero Mino time.")
    end
    for _ in 1:max_iter
        mid = 0.5 * (low + high)
        if _c1_lambda_of_r(c, mid) < target
            low = mid
        else
            high = mid
        end
    end
    return 0.5 * (low + high)
end

function _unsupported_capture(parameters, constants, outcome, reason)
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    return KerrGeoCapture(outcome.formula, outcome.energy_regime, outcome.outcome, parameters,
        constants, "Mino", (radial=outcome.roots, root_class=outcome.root_class),
        (kind=:future_horizon_regular_v, lambda_horizon=0.0),
        nothing, nothing,
        (radial=r -> kerr_radial_potential(a, energy, lz, q, r),
         polar_z=z -> kerr_polar_z_potential(a, energy, lz, q, z)),
        (radial=nothing, polar_z=nothing),
        (supported=false, reason=reason))
end

"""
    kerr_geo_capture(a, p, e, x; input=:apex, kwargs...)
    kerr_geo_capture(a, E, Lz, Q; input=:constants, polar_phase=0.0)
    kerr_geo_capture(a, constants; kwargs...)

Finite-window capture orbit (E ≥ 1, from infinity into the future horizon). Inclined
orbits need constants input; `polar_phase` is the polar argument at the horizon.
Constants without a capture component return a record with `Status.supported == false`.
APEX input whose periapsis p/(1 + e) lies outside the horizon throws an `ArgumentError`: a
capture orbit has no periapsis there.
"""
function kerr_geo_capture(a::Real, b::Real, c::Real, d::Real; input::Symbol=:apex,
        constants=nothing, polar_phase=nothing)
    if input === :constants
        constants_tuple = (E=b, Lz=c, Q=d)
        parameters = (a=a, input=:constants)
    elseif input === :apex
        # p/(1+e) is a root of R; outside the horizon it is a turning point that bounds every
        # motion from infinity (and below the separatrix the APEX constants may not exist)
        b / (1 + c) > _rplus(a) && throw(ArgumentError(
            "kerr_geo_capture: the periapsis p/(1+e) = $(b / (1 + c)) lies outside the horizon " *
            "r₊ = $(_rplus(a)); a capture orbit has no periapsis there. " *
            "Pass the constants of motion with input=:constants."))
        if constants === nothing
            c_ = kerr_geo_constants_of_motion(a, b, c, d)
            constants_tuple = (E=c_["E"], Lz=c_["Lz"], Q=c_["Q"])
        else
            constants_tuple = constants
        end
        parameters = (a=a, p=b, e=c, x=d, input=:apex)
    else
        error("Unknown capture input type. Use :apex or :constants.")
    end
    outcome = _outcome_at_infinity(a, constants_tuple.E, constants_tuple.Lz, constants_tuple.Q)
    outcome.outcome === :capture || return _unsupported_capture(parameters, constants_tuple,
        outcome, "Outcome is $(outcome.outcome), not capture.")
    return _capture_finite_window(parameters, constants_tuple, outcome;
        polar_phase=polar_phase === nothing ? 0.0 : polar_phase)
end

function kerr_geo_capture(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_capture(a, constants.E, constants.Lz, constants.Q; input=:constants, kwargs...)
end

function kerr_geo_capture(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_capture(a, constants[1], constants[2], constants[3]; input=:constants, kwargs...)
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoCapture)
    println(io, "KerrGeoCapture(")
    print(io, "    Formula = "); show(io, kg.Formula); println(io, ",")
    print(io, "    EnergyRegime = "); show(io, kg.EnergyRegime); println(io, ",")
    print(io, "    Outcome = "); show(io, kg.Outcome); println(io, ",")
    print(io, "    ConstantsOfMotion = "); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    ReferenceZero = "); show(io, kg.ReferenceZero); println(io, ",")
    print(io, "    Status = "); show(io, kg.Status); println(io, ",")
    print(io, ")")
end

# Radial closed forms (r(λ), λ(r)) of the two finite-window capture formulas:
# C3 (E > 1, two real roots inside the horizon plus a complex pair) and
# C1 (E = 1, one real root plus a complex pair).
function _capture_radial_model(formula, a, energy, lz, q)
    if formula === :hyperbolic_capture
        c = _c3_complex_parameters(a, energy, lz, q)
        c === nothing && return nothing
        return (params=c, rplus=c.rplus, lambda_of_r=r -> _c3_lambda_of_r(c, r),
            lambda_infinity=_c3_lambda_infinity(c),
            radius=λ -> _c3_radius_from_horizon_lambda(c, λ),
            rootdata=(radial=c.real_roots, complex=(rho=c.rho, eta=c.eta), lead=c.lead,
                      shape=_c3_shape(c)))
    else
        c = _c1_one_real_parameters(a, lz, q)
        c === nothing && return nothing
        return (params=c, rplus=c.rplus, lambda_of_r=r -> _c1_lambda_of_r(c, r),
            lambda_infinity=_c1_lambda_infinity(c),
            radius=λ -> _c1_radius_from_horizon_lambda(c, λ),
            rootdata=(radial=c.real_roots, complex=(rho=c.rho, eta=c.eta),
                      shape=_c1_shape(c)))
    end
end

"""
Finite-window capture orbit for E > 1 (two real roots inside the horizon plus a complex
pair) or E = 1 (one real root plus a complex pair): radial motion from infinity into the
future horizon, λ = 0 on the horizon where the regular v and ψ vanish. Equatorial orbits
and generic inclined ones with a fixed `polar_phase` are supported.
"""
function _capture_finite_window(parameters, constants, outcome; polar_phase=0.0)
    formula = outcome.formula
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    unsupported(msg) = _unsupported_capture(parameters, constants, outcome, msg)
    inclined = abs(q) > 1e-12
    inclined && parameters.input !== :constants && return unsupported(
        "Inclined capture needs constants input (the APEX polar-phase convention is not defined).")
    radial = _capture_radial_model(formula, a, energy, lz, q)
    radial === nothing && return unsupported(
        "The radial roots do not have the ordering required by the $formula formula.")
    polar = _window_polar_motion(a, energy, lz, q, polar_phase)

    lambda_horizon = radial.lambda_of_r(radial.rplus)
    lambda_infinity = lambda_horizon - radial.lambda_infinity
    R(r) = kerr_radial_potential(a, energy, lz, q, r)
    Θ(z) = kerr_polar_z_potential(a, energy, lz, q, z)
    r_of_lambda = radial.radius
    rdot(λ) = -sqrt(max(R(r_of_lambda(λ)), 0.0))
    # Mino time of radius r (λ = 0 on the horizon) and the polar increments between radii.
    function λ_of(r)
        r >= radial.rplus || throw(DomainError(r, "Capture radius must lie outside the future horizon."))
        return lambda_horizon - radial.lambda_of_r(r)
    end
    polar_phi(r1, r2) = polar.phi(λ_of(r1)) - polar.phi(λ_of(r2))
    polar_t(r1, r2) = polar.t(λ_of(r1)) - polar.t(λ_of(r2))
    # τ, v, ψ and the radial increments: radial spectral engine (infinity → horizon) +
    # polar primitive; τ, v, ψ vanish on the horizon (λ = 0)
    coords = _engine_coordinates(a, energy, lz, q, r_of_lambda, polar.primitive;
        domain=(lambda_infinity, 0.0), ends=(:infinity, :horizon), σ=-1.0,
        λ_bl=0.5 * lambda_infinity, λ_tau=0.0, λ_regular=0.0, σ_regular=-1.0)
    inc = _radius_increments(coords, λ_of, -1.0)
    radial_time, radial_phi, proper, radial_v, radial_psi = inc.radial_time_increment,
        inc.radial_phi_increment, inc.radial_proper_increment, inc.radial_v_increment,
        inc.radial_psi_increment
    tau, v, psi = coords.tau, coords.v, coords.psi
    rstar(λ) = λ >= -2e-13 ? NaN : kerr_rstar(a, r_of_lambda(λ))
    utheta(λ) = inclined ? -polar.uz(λ) / sqrt(max(1 - polar.z(λ)^2, 0.0)) : 0.0
    phase = _polar_phase_metadata(polar, polar_phase, :future_horizon_regular_endpoint)

    return KerrGeoCapture(formula, outcome.energy_regime, :capture, parameters, constants, "Mino",
        merge(radial.rootdata, (root_class=outcome.root_class,
                                polar=inclined ? polar.rootdata : nothing)),
        (kind=:future_horizon_regular_v, lambda_horizon=0.0,
         lambda_infinity=lambda_infinity, v_horizon=0.0, psi_horizon=0.0,
         rplus=radial.rplus, polar_phase=inclined ? polar_phase : nothing),
        (t=nothing, r=r_of_lambda, theta=polar.theta, phi=nothing, q=polar.z, tau=tau,
         rstar=rstar, u=nothing, v=v, psi=psi,
         radial_phi_increment=radial_phi, radial_time_increment=radial_time,
         radial_proper_increment=proper,
         total_phi_increment=(r1, r2) -> radial_phi(r1, r2) + polar_phi(r1, r2),
         total_time_increment=(r1, r2) -> radial_time(r1, r2) + polar_t(r1, r2),
         radial_psi_increment=radial_psi, radial_v_increment=radial_v,
         total_psi_increment=(r1, r2) -> radial_psi(r1, r2) + polar_phi(r1, r2),
         total_v_increment=(r1, r2) -> radial_v(r1, r2) + polar_t(r1, r2)),
        (ut=nothing, ur=rdot, utheta=utheta, uphi=nothing),
        (radial=R, polar_z=Θ),
        (radial=λ -> rdot(λ)^2 - R(r_of_lambda(λ)),
         polar_z=λ -> polar.uz(λ)^2 - Θ(polar.z(λ))),
        (supported=true, reason="ok", domain=_capture_domain(lambda_infinity),
         polar_phase=phase))
end

# ---------------------------------------------------------------------------------------
# Radial integrals over a C1 / C3 radial window by adaptive Gauss-Kronrod in
# u = log(r − r₊) (which absorbs the 1/(r − r₊) horizon pole of the BL rates). Used only
# by the closed-form members built outside the spectral engine: the E = 1 axis infall and
# the exact-extremal (a = 1) capture bases.
# ---------------------------------------------------------------------------------------
function _capture_sqrt_radial(c, r)
    complex_part = (r - c.rho)^2 + c.eta^2
    polynomial = haskey(c, :lead) ?
        c.lead * (r - c.r1) * (r - c.r2) * complex_part :
        2 * (r - c.x0) * complex_part
    return sqrt(max(polynomial, 0.0))
end

function _capture_radial_quadrature(integrand, c, r_left, r_right)
    left, right, sign = _capture_positive_radial_interval(r_left, r_right)
    left == right && return 0.0
    left >= c.rplus || error("Capture radial increment requires a window outside the future horizon.")
    rp = c.rplus
    function g(u)
        y = exp(u)
        r = rp + y
        return integrand(r) * y / _capture_sqrt_radial(c, r)
    end
    # left == r+ is allowed for integrands regular at the horizon (u -> -Inf)
    lower = left == rp ? -Inf : log(left - rp)
    # maxevals caps pathological inputs; a well-posed window converges in a few hundred
    value, _ = quadgk(g, lower, log(right - rp); rtol=1.0e-13, atol=0.0, maxevals=20_000)
    return sign * value
end

_capture_t_rate(a, energy, lz, r) =
    (r^2 + a^2) * kerr_radial_momentum(a, energy, lz, r) / kerr_delta(a, r)
_c1_radial_time_increment(a, lz, c, r_left, r_right) =
    _capture_radial_quadrature(r -> _capture_t_rate(a, 1.0, lz, r), c, r_left, r_right)
