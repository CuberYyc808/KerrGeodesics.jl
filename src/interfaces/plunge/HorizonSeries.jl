# Plunge reference API: near-horizon series of the ingoing/outgoing coordinates (series
# arithmetic in core/Series.jl) and the horizon anchoring of `kerr_geo_plunge`'s time origin.

function _rstar_series_constants(a)
    rp = _rplus(a)
    rm = _rminus(a)
    d = rp - rm
    alpha = 2 * rp / d
    cstar = rp - 2 * rp / d * log(oftype(rp, 2)) - 2 * rm / d * log(d / 2)
    return (rp=rp, rm=rm, d=d, alpha=alpha, cstar=cstar)
end

# The Taylor coefficients in x = r − r₊ of the polar part of dt/dλ, f = aLz − a²E(1 − z²),
# along the ingoing orbit as it reaches the horizon at λ_H: λ − λ_H = s(x) = −∫₀ˣ dx′/√R(r₊ + x′)
# from the series of R, z(λ_H + s) from z″ = Θ′(z)/2 = −(Q + Lz² + c) z + 2c z³ (c = a²(1 − E²))
# with z(λ_H) = z0 and dz/dλ(λ_H) = uz0, and f composed in series: no fit, every coefficient
# exact to rounding.
function _near_horizon_polar_series(a, energy, lz, q, z0, uz0, rseries, order)
    T = eltype(rseries)
    inverse_root = _series_inv(_series_sqrt_positive(rseries, order), order)
    s = zeros(T, order + 1)
    for n in 1:order
        s[n + 1] = -inverse_root[n] / n
    end
    c = -a^2 * _e2m1(energy)
    b = q + lz^2 + c
    zs = zeros(T, order + 1)
    zs[1] = z0
    order >= 1 && (zs[2] = uz0)
    for n in 0:order - 2
        cube = _series_mul(_series_mul(zs, zs, n), zs, n)
        zs[n + 3] = (-b * zs[n + 1] + 2c * cube[n + 1]) / ((n + 1) * (n + 2))
    end
    z = _series_compose(zs, s, order)
    f = a^2 * energy .* _series_mul(z, z, order)
    f[1] += a * lz - a^2 * energy
    return f
end

function _horizon_vq_coefficients(a, energy, lz, q, z0, uz0; order=10)
    order = max(1, order)
    T = _float_type(a, energy, lz, q)
    c = _rstar_series_constants(a)
    rp, rm, d, alpha = c.rp, c.rm, c.d, c.alpha

    a_series = zeros(T, order + 1)
    a_series[1] = 2 * rp
    order >= 1 && (a_series[2] = 2 * rp)
    order >= 2 && (a_series[3] = 1.0)

    delta = zeros(T, order + 1)
    order >= 1 && (delta[2] = d)
    order >= 2 && (delta[3] = 1.0)
    delta_regular = zeros(T, order + 1)
    delta_regular[1] = d
    order >= 1 && (delta_regular[2] = 1.0)

    pseries = zeros(T, order + 1)
    pseries[1] = energy * (rp^2 + a^2) - a * lz
    order >= 1 && (pseries[2] = 2 * energy * rp)
    order >= 2 && (pseries[3] = energy)
    pseries[1] <= 0 && error("The ingoing future-horizon series requires " *
        "P(r₊) = E(r₊² + a²) − a Lz > 0.")

    sseries = zeros(T, order + 1)
    sseries[1] = rp^2 + (a * energy - lz)^2 + q
    order >= 1 && (sseries[2] = 2 * rp)
    order >= 2 && (sseries[3] = 1.0)

    rseries = _series_mul(pseries, pseries, order) .- _series_mul(delta, sseries, order)
    fseries = _near_horizon_polar_series(a, energy, lz, q, z0, uz0, rseries, order)
    sqrt_r = _series_sqrt_positive(rseries, order)
    p_over_sqrt_r = _series_div(pseries, sqrt_r, order)
    one_minus = zeros(T, order + 1)
    one_minus[1] = 1.0
    one_minus .-= p_over_sqrt_r

    numerator = _series_mul(a_series, one_minus, order)
    shifted = zeros(T, order + 1)
    for n in 0:(order - 1)
        shifted[n + 1] = numerator[n + 2]
    end
    term1 = _series_div(shifted, delta_regular, order)
    term2 = -_series_div(fseries, sqrt_r, order)
    dvdr = term1 .+ term2

    vx = zeros(T, order + 1)
    for n in 1:order
        vx[n + 1] = dvdr[n] / n
    end

    rstar_power = zeros(T, order + 1)
    if order >= 1
        rstar_power[2] = (d^2 - 2 * rm) / d^2
    end
    for n in 2:order
        sign = iseven(n) ? 1 : -1
        rstar_power[n + 1] = sign * 2 * rm / (n * d^(n + 1))
    end
    qexp = _series_exp(rstar_power ./ alpha, order)
    q_of_x = zeros(T, order + 1)
    for n in 1:order
        q_of_x[n + 1] = qexp[n]
    end
    x_of_q = _series_revert_unit_linear(q_of_x, order)
    vq = _series_compose(vx, x_of_q, order)
    return vq[2:(order + 1)]
end

function _near_horizon_series_unavailable()
    nanfun(rstar; order=10) = NaN
    return (
        available=false,
        v_rstar=nanfun,
        u_rstar=nanfun,
        q_of_rstar=(rstar -> NaN),
        last_term_abs=(rstar; order=10) -> NaN,
        coefficients=Float64[],
        vH=NaN,
        metadata=(
            rstar_convention="kerr_rstar",
            near_horizon_branch="none",
            P_plus_sign="not_evaluated",
            near_horizon_series_order=0,
            near_horizon_series_variable="q_exp_rstar_minus_Cstar_over_alpha",
            lambda_to_rstar_path="kerr_rstar_of_trajectory_radius",
        ),
    )
end

function _build_near_horizon_uv_series(a, energy, lz, q, theta, uz, v, lambda_of_radius,
        lambda_r0; order=10, horizon_offset=1e-4)
    order = max(1, order)
    c = _rstar_series_constants(a)
    rp, alpha, cstar = c.rp, c.alpha, c.cstar
    pplus = energy * (rp^2 + a^2) - a * lz
    if pplus <= 0
        return _near_horizon_series_unavailable()
    end

    # the polar position and velocity where the orbit reaches the horizon
    lambda_h = lambda_of_radius(rp) - lambda_r0
    coeffs = _horizon_vq_coefficients(a, energy, lz, q, cos(theta(lambda_h)), uz(lambda_h);
        order=order)

    anchors = typeof(rp)[]
    for scale in (1.0, 1.25, 1.5, 2.0, 3.0)
        x = horizon_offset * scale
        radius = rp + x
        lambda = lambda_of_radius(radius) - lambda_r0
        rstar = kerr_rstar(a, radius)
        qvar = exp((rstar - cstar) / alpha)
        correction = zero(rp)
        for n in 1:order
            correction += coeffs[n] * qvar^n
        end
        push!(anchors, v(lambda) - correction)
    end
    vH = sum(anchors) / length(anchors)

    q_of_rstar(rstar) = exp((rstar - cstar) / alpha)
    function v_rstar(rstar; order=length(coeffs))
        nmax = max(1, min(order, length(coeffs)))
        qvar = q_of_rstar(rstar)
        value = vH
        for n in 1:nmax
            value += coeffs[n] * qvar^n
        end
        return value
    end
    u_rstar(rstar; order=length(coeffs)) = v_rstar(rstar; order=order) - 2 * rstar
    function last_term_abs(rstar; order=length(coeffs))
        nmax = max(1, min(order, length(coeffs)))
        return abs(coeffs[nmax] * q_of_rstar(rstar)^nmax)
    end

    return (
        available=true,
        v_rstar=v_rstar,
        u_rstar=u_rstar,
        q_of_rstar=q_of_rstar,
        last_term_abs=last_term_abs,
        coefficients=coeffs,
        vH=vH,
        metadata=(
            rstar_convention="kerr_rstar",
            near_horizon_branch="ingoing_future_horizon",
            P_plus_sign="positive_required",
            near_horizon_series_order=order,
            near_horizon_series_variable="q_exp_rstar_minus_Cstar_over_alpha",
            lambda_to_rstar_path="kerr_rstar_of_trajectory_radius",
        ),
    )
end

function _estimate_direct_future_horizon_v_anchor(
        a, t, lambda_of_radius, lambda_r0;
        horizon_offset=1e-4,
        scales=(1.0, 1.25, 1.5, 2.0, 3.0))
    rp = _rplus(a)
    anchors = typeof(rp)[]
    for scale in scales
        radius = rp + horizon_offset * scale
        lambda = lambda_of_radius(radius) - lambda_r0
        push!(anchors, t(lambda) + kerr_rstar(a, radius))
    end
    return sum(anchors) / length(anchors), anchors
end

function _estimate_future_horizon_v_anchor(
        a, energy, lz, q, root_class, theta, uz, t, v, lambda_of_radius, lambda_r0;
        order=10,
        horizon_offset=1e-4)
    if root_class == "Real2"
        series = _build_near_horizon_uv_series(
            a, energy, lz, q, theta, uz, v, lambda_of_radius, lambda_r0;
            order=order,
            horizon_offset=horizon_offset,
        )
        if series.available && isfinite(series.vH)
            return (
                vH=series.vH,
                method="real2_horizon_series_overlap_anchor",
                direct_anchor_values=Float64[],
            )
        end
    end
    direct_vH, anchors = _estimate_direct_future_horizon_v_anchor(
        a, t, lambda_of_radius, lambda_r0;
        horizon_offset=horizon_offset,
    )
    return (
        vH=direct_vH,
        method="direct_t_plus_rstar_overlap_anchor",
        direct_anchor_values=anchors,
    )
end
