# Elementary pieces shared by the closed-form radial models (Plunge, Capture, Critical, axis
# infall, extremal): the quadratic-denominator primitives and the Legendre J2 combination, plus
# the Mino-time finiteness check and the real-root radii used by the member tracks.

_finite_mino(λ) = (isfinite(λ) || throw(DomainError(λ, "Mino time must be finite.")); float(λ))

_root_radii(classification) = Tuple(
    item.radius for item in classification.Status.root_structure.real_roots)

function _real_atanh(x)
    abs(x) == 1 && return copysign(Inf, x)
    return abs(x) < 1 ? atanh(x) : atanh(inv(x))
end

function _j_inv(k, y)
    if iszero(k)
        return y
    elseif k > 0
        return _real_atanh(sqrt(k) * y) / sqrt(k)
    end
    return atan(sqrt(-k) * y) / sqrt(-k)
end

# ∫ dy / (a0 + b0 y²)
function _quadratic_denominator_primitive(a0, b0, y)
    b0 > 0 || error("The quadratic denominator requires a positive quadratic coefficient.")
    iszero(a0) && return -inv(b0 * y)
    if a0 > 0
        return atan(y * sqrt(b0 / a0)) / sqrt(a0 * b0)
    end
    y0 = sqrt(-a0 / b0)
    return -_real_atanh(y / y0) / (b0 * y0)
end

