# Truncated power-series arithmetic (product, inverse, quotient, square root, exponential,
# composition and reversion) used by the near-horizon series of the coordinates, in the
# floating-point type of the coefficients.

_series_type(a) = float(eltype(a))

function _series_mul(a, b, order)
    c = zeros(promote_type(_series_type(a), _series_type(b)), order + 1)
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
    T = _series_type(a)
    b = zeros(T, order + 1)
    b[1] = inv(a[1])
    for n in 1:order
        s = zero(T)
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
    T = _series_type(a)
    b = zeros(T, order + 1)
    b[1] = sqrt(a[1])
    for n in 1:order
        s = zero(T)
        for k in 1:(n - 1)
            s += b[k + 1] * b[n - k + 1]
        end
        b[n + 1] = (a[n + 1] - s) / (2 * b[1])
    end
    return b
end

function _series_exp(a, order)
    T = _series_type(a)
    b = zeros(T, order + 1)
    b[1] = exp(a[1])
    for n in 1:order
        s = zero(T)
        for k in 1:n
            s += k * a[k + 1] * b[n - k + 1]
        end
        b[n + 1] = s / n
    end
    return b
end

function _series_compose(f, g, order)
    T = promote_type(_series_type(f), _series_type(g))
    if abs(g[1]) > 100 * eps(T)
        error("Series composition expects an inner series with zero constant term.")
    end
    out = zeros(T, order + 1)
    power = zeros(T, order + 1)
    power[1] = 1
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
    T = _series_type(f)
    if abs(f[1]) > 100 * eps(T) || abs(f[2] - 1) > _tol(T, 1e-10)
        error("Series reversion requires q = x + O(x²).")
    end
    g = zeros(T, order + 1)
    g[2] = 1
    for n in 2:order
        trial = copy(g)
        trial[n + 1] = 0
        composed = _series_compose(f, trial, order)
        g[n + 1] = -composed[n + 1]
    end
    return g
end
