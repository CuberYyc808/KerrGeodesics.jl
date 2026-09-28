# Mino-time four-velocity, dτ/dλ and the geodesic residuals of a member, shared by the member
# constructors: each member supplies only its r(λ), polar position z(λ), dz/dλ and radial
# direction.

"""
    _kinematics(a, E, Lz, Q; r, rbl=r, z, uz, sin2, sign_r, R)

Closures of λ for one member. `r` is the radius on the member's whole domain, `rbl` the
radius restricted to where Boyer-Lindquist t, φ exist (it throws outside), `z`, `uz` the polar
position and velocity, `sin2` = sin²θ from the polar solution (formed without cancellation:
1 − z² from a rounded z loses the Lz/sin²θ rate of near-axis orbits), `sign_r(λ)` the sign of
dr/dλ, `R(r)` the radial potential. The t and φ rates are the radial engine's
(`_plain_rates`); the polar φ rate vanishes identically for Lz = 0.
"""
function _kinematics(a, E, L, Q; r, rbl=r, z, uz, sin2, sign_r, R)
    c = _rc(a, E, L, Q)
    Θ(zv) = kerr_polar_z_potential(a, E, L, Q, zv)
    ur(λ) = sign_r(λ) * sqrt(max(R(r(λ)), 0.0))
    function utheta(λ)
        s = sqrt(sin2(λ))
        return iszero(s) ? NaN : -uz(λ) / s            # (on the axis θ is a chart pole)
    end
    ut(λ) = _plain_rates(c, rbl(λ))[1] + a * L - a^2 * E * sin2(λ)
    uphi(λ) = _plain_rates(c, rbl(λ))[2] + (iszero(L) ? 0.0 : L / sin2(λ))
    dtau_dlambda(λ) = r(λ)^2 + a^2 * z(λ)^2
    # g_{μν} u^μ u^ν + 1 with the proper-time velocity u = (dx/dλ)/Σ
    function normalization_residual(λ)
        rv = rbl(λ); zv = z(λ)
        Σ = rv^2 + a^2 * zv^2
        Δ = kerr_delta(a, rv)
        s2 = sin2(λ)
        vt, vr, vθ, vφ = ut(λ) / Σ, ur(λ) / Σ, utheta(λ) / Σ, uphi(λ) / Σ
        return -(1 - 2rv / Σ) * vt^2 - 2 * (2a * rv * s2 / Σ) * vt * vφ + Σ / Δ * vr^2 +
            Σ * vθ^2 + s2 * (rv^2 + a^2 + 2a^2 * rv * s2 / Σ) * vφ^2 + 1
    end
    return (ur=ur, uz=uz, utheta=utheta, ut=ut, uphi=uphi, dtau_dlambda=dtau_dlambda,
        radial_residual=λ -> ur(λ)^2 - R(r(λ)), polar_residual=λ -> uz(λ)^2 - Θ(z(λ)),
        normalization_residual=normalization_residual, R=R, Θ=Θ)
end
