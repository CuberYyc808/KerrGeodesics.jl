# Mino-time four-velocity, dτ/dλ and the geodesic residuals of a member, shared by the member
# constructors: each member supplies only its r(λ), polar position z(λ), dz/dλ and radial
# direction.

"""
    _kinematics(a, E, Lz, Q, r, rbl, position, sign_r, R)

Closures of λ for one member. `r` is the radius on the member's whole domain, `rbl` the
radius restricted to where Boyer-Lindquist t, φ exist (it throws outside), `position(λ)` the
polar `(z, dz/dλ, sin²θ)`, sin²θ from the polar solution (formed without cancellation:
1 − z² from a rounded z loses the Lz/sin²θ rate of near-axis orbits), `sign_r(λ)` the sign of
dr/dλ, `R(r)` the radial potential. The t and φ rates are the radial engine's
(`_plain_rates`); the polar φ rate vanishes identically for Lz = 0. Every closure captures one
`_KinematicState` and calls a named function on it (`_kin_ut`, …), so no closure's type
contains another's.
"""
function _kinematics(a, E, L, Q, r, rbl, position, sign_r, R; radial_track=nothing)
    k = _KinematicState(a, E, L, Q, _rc(a, E, L, Q, R), r, rbl, position, sign_r, R,radial_track)
    Θ(zv) = kerr_polar_z_potential(a, E, L, Q, zv)
    return (ur=λ -> _kin_ur(k, λ), uz=λ -> _kin_uz(k, λ), utheta=λ -> _kin_utheta(k, λ),
        ut=λ -> _kin_ut(k, λ), uphi=λ -> _kin_uphi(k, λ), dtau_dlambda=λ -> _kin_dtau(k, λ),
        radial_residual=λ -> _kin_radial_residual(k,λ),
        polar_residual=λ -> _kin_uz(k, λ)^2 - Θ(_kin_z(k, λ)),
        normalization_residual=λ -> _kin_normalization(k, λ), R=R, Θ=Θ, state=k)
end

struct _KinematicState{T,C,Fr,Fb,Fp,Fg,FR,Ft}
    a::T; E::T; L::T; Q::T
    c::C                                     # _RadialConstants for the t, φ rates
    r::Fr; rbl::Fb; position::Fp; sign_r::Fg; R::FR
    radial_track::Ft
end
_KinematicState(a, E, L, Q, c, r, rbl, position, sign_r, R,radial_track) =
    _KinematicState(promote(a, E, L, Q)..., c, r, rbl, position, sign_r, R,radial_track)

_kin_z(k::_KinematicState, λ) = k.position(λ)[1]
_kin_uz(k::_KinematicState, λ) = k.position(λ)[2]
_kin_sin2(k::_KinematicState, λ) = k.position(λ)[3]

function _kin_ur(k::_KinematicState,λ)
    state=_radial_state(k.radial_track,λ)
    return state === nothing ? k.sign_r(λ)*sqrt(max(k.R(k.r(λ)),0.0)) : state.velocity
end
function _kin_radial_residual(k::_KinematicState,λ)
    state=_radial_state(k.radial_track,λ)
    potential=state === nothing ? k.R(k.r(λ)) : _wide_evalpoly(state.gap,state.chart.shifted)
    return _kin_ur(k,λ)^2-potential
end
function _kin_radial_rates(k::_KinematicState,λ)
    radius=k.rbl(λ)
    state=_radial_state(k.radial_track,λ)
    return state === nothing ? _plain_rates(k.c,radius) : _relative_rates(k.c,state,:plain,1.0)
end
function _kin_utheta(k::_KinematicState, λ)
    s = sqrt(_kin_sin2(k, λ))
    return iszero(s) ? NaN : -_kin_uz(k, λ) / s            # (on the axis θ is a chart pole)
end
_kin_ut(k::_KinematicState, λ) = _kin_radial_rates(k,λ)[1] + k.a * k.L - k.a^2 * k.E * _kin_sin2(k, λ)
_kin_uphi(k::_KinematicState, λ) = _kin_radial_rates(k,λ)[2] + (iszero(k.L) ? 0.0 : k.L / _kin_sin2(k, λ))
_kin_dtau(k::_KinematicState, λ) = k.r(λ)^2 + k.a^2 * _kin_z(k, λ)^2
# g_{μν} u^μ u^ν + 1 with the proper-time velocity u = (dx/dλ)/Σ
function _kin_normalization(k::_KinematicState, λ)
    a = k.a
    rv = k.rbl(λ); zv = _kin_z(k, λ)
    Σ = rv^2 + a^2 * zv^2
    s2 = _kin_sin2(k, λ)
    state=_radial_state(k.radial_track,λ)
    if state !== nothing
        # Substitution of the separated rates into g(u,u)+1 leaves the radial
        # and polar equation residuals; this avoids cancelling divergent BL terms.
        delta=state.gap*(state.gap+state.chart.separation)
        polar=_kin_uz(k,λ)^2-kerr_polar_z_potential(k.a,k.E,k.L,k.Q,zv)
        return (_kin_radial_residual(k,λ)/delta+polar/s2)/Σ
    end
    Δ = kerr_delta(a, rv)
    vt, vr, vθ, vφ = _kin_ut(k, λ) / Σ, _kin_ur(k, λ) / Σ, _kin_utheta(k, λ) / Σ, _kin_uphi(k, λ) / Σ
    return -(1 - 2rv / Σ) * vt^2 - 2 * (2a * rv * s2 / Σ) * vt * vφ + Σ / Δ * vr^2 +
        Σ * vθ^2 + s2 * (rv^2 + a^2 + 2a^2 * rv * s2 / Σ) * vφ^2 + 1
end
