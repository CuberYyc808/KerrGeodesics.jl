module InitialConditions

include("OrbitClass.jl")
using .OrbitClass
include("OrbitalDuration.jl")
using .OrbitalDuration
using Elliptic

export _theta_phase_from_theta, _radial_phase_from_options

function _theta_phase_from_theta(a, energy, lz, q, theta)
    zm, zp = polar_roots(a, energy, lz, q)
    if zm <= 0 || iszero(a) || iszero(q)
        return 0.0
    end
    z = cos(theta)
    amplitude = z / sqrt(zm)
    if abs(amplitude) > 1 + 1e-12
        error("initial_theta is outside the allowed polar range for this plunge.")
    end
    amplitude = clamp(amplitude, -1.0, 1.0)
    xi_theta = sqrt(a^2 * (1 - energy^2) * zp)
    ktheta = zm / zp
    return Elliptic.F(asin(amplitude), ktheta) / xi_theta
end

function _radial_phase_from_options(a, energy, lz, q; radial_phase=nothing, initial_radius=nothing, radial_start=:outer_turning)
    if radial_phase !== nothing
        return radial_phase
    end
    lambda_max, _, lambda_of_radius = lambda_of_r(a, energy, lz, q)
    if initial_radius !== nothing
        phase = lambda_of_radius(initial_radius)
        if phase === nothing
            error("initial_radius is outside the implemented plunge radial range.")
        end
        return phase
    end
    if radial_start == :outer_turning
        return 0.0
    elseif radial_start == :inner_turning
        return lambda_max
    else
        error("Unknown radial_start. Use :outer_turning or :inner_turning.")
    end
end

end
