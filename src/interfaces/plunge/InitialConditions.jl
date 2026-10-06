# Plunge reference API: radial and polar phases from the user's initial data.

function _theta_phase_from_theta(a, energy, lz, q, theta)
    zm, xi_theta, ktheta = _plunge_polar_parameters(a, energy, lz, q)
    if zm <= 0 || iszero(q)
        return zero(zm)
    end
    z = cos(theta)
    amplitude = z / sqrt(zm)
    if abs(amplitude) > 1 + _tol(typeof(amplitude), 1e-12)
        error("initial_theta is outside the allowed polar range for this plunge.")
    end
    amplitude = clamp(amplitude, -1.0, 1.0)
    return _F(asin(amplitude), ktheta) / xi_theta
end

function _radial_phase_from_options(a, energy, lz, q; radial_phase=nothing, initial_radius=nothing, radial_start=:turning_point)
    if radial_phase !== nothing
        return radial_phase
    end
    lambda_max, _, lambda_of_radius = lambda_of_r(a, energy, lz, q)
    if initial_radius !== nothing
        phase = lambda_of_radius(initial_radius)
        if phase === nothing
            error("initial_radius = $initial_radius is outside the radial range of this plunge.")
        end
        return phase
    end
    if radial_start in (:turning_point, :outer_turning)
        return zero(_float_type(a, energy, lz, q))
    elseif radial_start == :inner_turning
        return lambda_max
    else
        error("Unknown radial_start = $(repr(radial_start)). Use :turning_point (or " *
            ":outer_turning) for the outer turning point, or :inner_turning for the inner end " *
            "of the radial range.")
    end
end
