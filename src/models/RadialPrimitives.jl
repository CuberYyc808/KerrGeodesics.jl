# Elementary and elliptic pieces shared by the closed-form radial models (Plunge, Capture,
# Critical, axis infall, extremal), plus the Mino-time finiteness check and the real-root
# radii used by the member tracks.

_finite_mino(λ) = (isfinite(λ) || throw(DomainError(λ, "Mino time must be finite.")); float(λ))

_root_radii(classification) = Tuple(
    item.radius for item in classification.Status.root_structure.real_roots)

function _real_atanh(x)
    abs(x) == 1 && return copysign(Inf, x)
    return 0.5 * log(abs((1 + x) / (1 - x)))
end

function _j_inv(k, y)
    if abs(k) <= 1.0e-15
        return y
    elseif k > 0
        return _real_atanh(sqrt(k) * y) / sqrt(k)
    end
    return atan(sqrt(-k) * y) / sqrt(-k)
end

# ∫ dy / (a0 + b0 y²)
function _quadratic_denominator_primitive(a0, b0, y)
    scale = max(1.0, abs(a0), abs(b0))
    abs(a0) > 16 * eps(Float64) * scale || error(
        "The pole coincides with a radial root (a0 = 0 in ∫ dy/(a0 + b0 y²)).")
    b0 > 0 || error("The quadratic denominator requires a positive quadratic coefficient.")
    if a0 > 0
        return atan(y * sqrt(b0 / a0)) / sqrt(a0 * b0)
    end
    y0 = sqrt(-a0 / b0)
    return log(abs((y - y0) / (y + y0))) / (2 * b0 * y0)
end

function _pi_real(n, phi, m)
    abs(n - 1) <= 8 * eps(Float64) && error(
        "Elliptic Pi characteristic is at its separate n=1 limit.")
    return n > 1 ? elliptic_pi(n, phi, m) : Elliptic.Pi(n, phi, m)
end

function _j2_legendre(n, m, phi)
    boundary = sin(phi) * cos(phi) *
        sqrt(max(1 - m * sin(phi)^2, 0.0)) /
        (1 - n * sin(phi)^2)
    acoef = 1 / (2 * (n - 1))
    bcoef = n / (2 * (m - n) * (n - 1))
    ccoef = (2 * m * n - 3 * m - n^2 + 2 * n) /
        (2 * (m - n) * (n - 1))
    dcoef = -n^2 / (2 * (m - n) * (n - 1))
    return acoef * Elliptic.F(phi, m) +
           bcoef * Elliptic.E(phi, m) +
           ccoef * _pi_real(n, phi, m) + dcoef * boundary
end
