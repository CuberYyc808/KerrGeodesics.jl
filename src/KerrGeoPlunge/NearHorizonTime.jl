module NearHorizonTime

export kerr_rstar,
       _build_near_horizon_uv_series,
       _near_horizon_series_unavailable,
       _estimate_future_horizon_v_anchor

function _rplus(a)
    return 1 + sqrt(1 - a^2)
end

function _rminus(a)
    return 1 - sqrt(1 - a^2)
end

"""
    kerr_rstar(a, r)

Return the Kerr tortoise coordinate in the package convention for exterior
radii `r > rplus`. The function returns `NaN` for radii at or inside the outer
horizon.
"""
function kerr_rstar(a::Real, r::Real)
    rp = _rplus(a)
    rm = _rminus(a)
    r <= rp && return NaN
    return r + 2 * rp / (rp - rm) * log((r - rp) / 2) -
           2 * rm / (rp - rm) * log((r - rm) / 2)
end

function _rstar_series_constants(a)
    rp = _rplus(a)
    rm = _rminus(a)
    d = rp - rm
    alpha = 2 * rp / d
    cstar = rp - 2 * rp / d * log(2) - 2 * rm / d * log(d / 2)
    return (rp=rp, rm=rm, d=d, alpha=alpha, cstar=cstar)
end

function _series_mul(a, b, order)
    c = zeros(Float64, order + 1)
    for i in 0:order
        ai = a[i + 1]
        iszero(ai) && continue
        for j in 0:(order - i)
            c[i + j + 1] += ai * b[j + 1]
        end
    end
    return c
end

function _series_inv(a, order)
    if iszero(a[1])
        error("Series inverse requires a nonzero constant term.")
    end
    b = zeros(Float64, order + 1)
    b[1] = inv(a[1])
    for n in 1:order
        s = 0.0
        for k in 1:n
            s += a[k + 1] * b[n - k + 1]
        end
        b[n + 1] = -s / a[1]
    end
    return b
end

function _series_div(a, b, order)
    return _series_mul(a, _series_inv(b, order), order)
end

function _series_sqrt_positive(a, order)
    if a[1] <= 0
        error("Positive-root series square root requires positive leading coefficient.")
    end
    b = zeros(Float64, order + 1)
    b[1] = sqrt(a[1])
    for n in 1:order
        s = 0.0
        for k in 1:(n - 1)
            s += b[k + 1] * b[n - k + 1]
        end
        b[n + 1] = (a[n + 1] - s) / (2 * b[1])
    end
    return b
end

function _series_exp(a, order)
    b = zeros(Float64, order + 1)
    b[1] = exp(a[1])
    for n in 1:order
        s = 0.0
        for k in 1:n
            s += k * a[k + 1] * b[n - k + 1]
        end
        b[n + 1] = s / n
    end
    return b
end

function _series_compose(f, g, order)
    if abs(g[1]) > 100 * eps(Float64)
        error("Series composition expects an inner series with zero constant term.")
    end
    out = zeros(Float64, order + 1)
    power = zeros(Float64, order + 1)
    power[1] = 1.0
    for n in 0:order
        if n > 0
            power = _series_mul(power, g, order)
        end
        if !iszero(f[n + 1])
            out .+= f[n + 1] .* power
        end
    end
    return out
end

function _series_revert_unit_linear(f, order)
    if abs(f[1]) > 100 * eps(Float64) || abs(f[2] - 1) > 1e-10
        error("B50 q(x) reversion expects q=x+O(x^2).")
    end
    g = zeros(Float64, order + 1)
    g[2] = 1.0
    for n in 2:order
        trial = copy(g)
        trial[n + 1] = 0.0
        composed = _series_compose(f, trial, order)
        g[n + 1] = -composed[n + 1]
    end
    return g
end

function _b50_vq_coefficients(a, energy, lz, q, fcoeffs; order=10)
    order = max(1, min(order, 10))
    c = _rstar_series_constants(a)
    rp, rm, d, alpha = c.rp, c.rm, c.d, c.alpha

    a_series = zeros(Float64, order + 1)
    a_series[1] = 2 * rp
    order >= 1 && (a_series[2] = 2 * rp)
    order >= 2 && (a_series[3] = 1.0)

    delta = zeros(Float64, order + 1)
    order >= 1 && (delta[2] = d)
    order >= 2 && (delta[3] = 1.0)
    delta_regular = zeros(Float64, order + 1)
    delta_regular[1] = d
    order >= 1 && (delta_regular[2] = 1.0)

    pseries = zeros(Float64, order + 1)
    pseries[1] = energy * (rp^2 + a^2) - a * lz
    order >= 1 && (pseries[2] = 2 * energy * rp)
    order >= 2 && (pseries[3] = energy)
    pseries[1] <= 0 && error("B50 ingoing future-horizon series requires P_plus > 0.")

    sseries = zeros(Float64, order + 1)
    sseries[1] = rp^2 + (a * energy - lz)^2 + q
    order >= 1 && (sseries[2] = 2 * rp)
    order >= 2 && (sseries[3] = 1.0)

    fseries = zeros(Float64, order + 1)
    for n in 0:min(order, length(fcoeffs) - 1)
        fseries[n + 1] = fcoeffs[n + 1]
    end

    rseries = _series_mul(pseries, pseries, order) .- _series_mul(delta, sseries, order)
    sqrt_r = _series_sqrt_positive(rseries, order)
    p_over_sqrt_r = _series_div(pseries, sqrt_r, order)
    one_minus = zeros(Float64, order + 1)
    one_minus[1] = 1.0
    one_minus .-= p_over_sqrt_r

    numerator = _series_mul(a_series, one_minus, order)
    shifted = zeros(Float64, order + 1)
    for n in 0:(order - 1)
        shifted[n + 1] = numerator[n + 2]
    end
    term1 = _series_div(shifted, delta_regular, order)
    term2 = -_series_div(fseries, sqrt_r, order)
    dvdr = term1 .+ term2

    vx = zeros(Float64, order + 1)
    for n in 1:order
        vx[n + 1] = dvdr[n] / n
    end

    rstar_power = zeros(Float64, order + 1)
    if order >= 1
        rstar_power[2] = (d^2 - 2 * rm) / d^2
    end
    for n in 2:order
        sign = iseven(n) ? 1.0 : -1.0
        rstar_power[n + 1] = sign * 2 * rm / (n * d^(n + 1))
    end
    qexp = _series_exp(rstar_power ./ alpha, order)
    q_of_x = zeros(Float64, order + 1)
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
            rstar_convention="code_log_halves",
            near_horizon_branch="unavailable_for_current_branch",
            P_plus_sign="not_evaluated",
            near_horizon_series_order=0,
            near_horizon_series_variable="q_exp_rstar_minus_Cstar_over_alpha",
            u_status="direct_trajectory_callable_only",
            v_status="direct_trajectory_callable_only",
            lambda_to_rstar_path="numerical_or_existing_trajectory_not_replaced",
            series_switch_tolerance_target=1e-15,
            series_switch_status="not_configured",
            regular_time_coefficients_status="not_configured",
        ),
    )
end

function _solve_dense_linear_system(matrix, rhs)
    n = length(rhs)
    a = [matrix[i, j] for i in 1:n, j in 1:n]
    b = copy(rhs)
    for k in 1:n
        pivot = k
        pivot_abs = abs(a[k, k])
        for i in (k + 1):n
            candidate = abs(a[i, k])
            if candidate > pivot_abs
                pivot = i
                pivot_abs = candidate
            end
        end
        pivot_abs <= eps(Float64) && error("Singular near-horizon regular-time fit system.")
        if pivot != k
            for j in k:n
                a[k, j], a[pivot, j] = a[pivot, j], a[k, j]
            end
            b[k], b[pivot] = b[pivot], b[k]
        end
        for i in (k + 1):n
            factor = a[i, k] / a[k, k]
            a[i, k] = 0.0
            for j in (k + 1):n
                a[i, j] -= factor * a[k, j]
            end
            b[i] -= factor * b[k]
        end
    end
    x = zeros(Float64, n)
    for k in n:-1:1
        s = b[k]
        for j in (k + 1):n
            s -= a[k, j] * x[j]
        end
        x[k] = s / a[k, k]
    end
    return x
end

function _fit_near_horizon_regular_time_coefficients(
        a, energy, lz, theta, lambda_of_radius, lambda_r0, rp;
        order=10,
        horizon_offset=1e-4)
    fit_order = min(order, 10)
    scales = [1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 15.0]
    scales = scales[1:(fit_order + 1)]
    vandermonde = zeros(Float64, length(scales), fit_order + 1)
    values = zeros(Float64, length(scales))
    for (i, scale) in pairs(scales)
        radius = rp + horizon_offset * scale
        lambda = lambda_of_radius(radius) - lambda_r0
        z = cos(theta(lambda))
        values[i] = -a^2 * energy * (1 - z^2) + a * lz
        for j in 0:fit_order
            vandermonde[i, j + 1] = scale^j
        end
    end
    scaled_coeffs = _solve_dense_linear_system(vandermonde, values)
    fcoeffs = zeros(Float64, order + 1)
    for j in 0:fit_order
        fcoeffs[j + 1] = scaled_coeffs[j + 1] / horizon_offset^j
    end
    return fcoeffs, fit_order
end

function _build_near_horizon_uv_series(a, energy, lz, q, theta, t, lambda_of_radius, lambda_r0;
        order=10,
        horizon_offset=1e-4)
    order = max(1, min(order, 10))
    c = _rstar_series_constants(a)
    rp, alpha, cstar = c.rp, c.alpha, c.cstar
    pplus = energy * (rp^2 + a^2) - a * lz
    if pplus <= 0
        return _near_horizon_series_unavailable()
    end

    fcoeffs, fit_order = _fit_near_horizon_regular_time_coefficients(
        a, energy, lz, theta, lambda_of_radius, lambda_r0, rp;
        order=order,
        horizon_offset=horizon_offset,
    )
    coeffs = _b50_vq_coefficients(a, energy, lz, q, fcoeffs; order=order)

    anchors = Float64[]
    for scale in (1.0, 1.25, 1.5, 2.0, 3.0)
        x = horizon_offset * scale
        radius = rp + x
        lambda = lambda_of_radius(radius) - lambda_r0
        rstar = kerr_rstar(a, radius)
        qvar = exp((rstar - cstar) / alpha)
        correction = 0.0
        for n in 1:order
            correction += coeffs[n] * qvar^n
        end
        push!(anchors, t(lambda) + rstar - correction)
    end
    vH = sum(anchors) / length(anchors)

    q_of_rstar(rstar) = exp((rstar - cstar) / alpha)
    function v_rstar(rstar; order=10)
        nmax = max(1, min(order, length(coeffs)))
        qvar = q_of_rstar(rstar)
        value = vH
        for n in 1:nmax
            value += coeffs[n] * qvar^n
        end
        return value
    end
    u_rstar(rstar; order=10) = v_rstar(rstar; order=order) - 2 * rstar
    function last_term_abs(rstar; order=10)
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
            rstar_convention="code_log_halves",
            near_horizon_branch="ingoing_future_horizon",
            P_plus_sign="positive_required",
            near_horizon_series_order=order,
            near_horizon_series_variable="q_exp_rstar_minus_Cstar_over_alpha",
            u_status="diverges_linearly_for_ingoing_future_horizon",
            v_status="finite_regular_horizon_coordinate",
            lambda_to_rstar_path="numerical_or_existing_trajectory_not_replaced",
            series_switch_tolerance_target=1e-15,
            series_switch_status="diagnostic_callable_available_no_silent_production_switch",
            regular_time_coefficients_status="local_cutoff_F0_to_F$(fit_order)_interpolation_from_existing_theta_trajectory_higher_Fn_zero",
        ),
    )
end

function _estimate_direct_future_horizon_v_anchor(
        a, t, lambda_of_radius, lambda_r0;
        horizon_offset=1e-4,
        scales=(1.0, 1.25, 1.5, 2.0, 3.0))
    rp = 1 + sqrt(1 - a^2)
    anchors = Float64[]
    for scale in scales
        radius = rp + horizon_offset * scale
        lambda = lambda_of_radius(radius) - lambda_r0
        push!(anchors, t(lambda) + kerr_rstar(a, radius))
    end
    return sum(anchors) / length(anchors), anchors
end

function _estimate_future_horizon_v_anchor(
        a, energy, lz, q, root_class, theta, t, lambda_of_radius, lambda_r0;
        order=10,
        horizon_offset=1e-4)
    if root_class == "Real2"
        series = _build_near_horizon_uv_series(
            a, energy, lz, q, theta, t, lambda_of_radius, lambda_r0;
            order=order,
            horizon_offset=horizon_offset,
        )
        if series.available && isfinite(series.vH)
            return (
                vH=series.vH,
                method="real2_B50_series_overlap_anchor",
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

end
