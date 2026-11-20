# APEX reference API: the orbit of `kerr_geo_orbit` as a typed record evaluated by its radial
# and polar phases.
#
# Every oscillating part is written with the Jacobi functions (s, c, d) = (sn, cn, dn)(y) of
# the phase reduced to |y| ≤ K. With F(am y | m) = y the Legendre forms lose their secular
# parts exactly,
#
#     Π(n|m) y/K − Π(n; am y|m) = n (P(n|m) y/K − P(n; am y|m)),
#     E(m) y/K − E(am y|m)      = −m (D(m) y/K − D(am y|m)),
#
# with the pole part P = (Π − F)/n = R_J/3 and D = (F − E)/m = R_D/3. Each increment is then
# a multiple of n or m, so it vanishes exactly where they do (e = 0, |x| = 1), and the
# complement 1 − n is passed as it is formed without cancellation (x² for the polar
# characteristic 1 − x²). The phase is reduced at the turning points, q = jπ, so that the
# turning points themselves are evaluated exactly.

# One radial or polar phase: (sn, cn, dn) of the reduced argument y, |y| ≤ K, whether an odd
# number of half periods 2K was removed (sn and cn change sign), and y. The phase is offset from
# its nearest turning point, q = jπ, which lies at u = jK (radial, u = Kq/π) or at u = (2j + 1)K
# (polar, u = K(2q/π + 1)); `scale` maps the offset to v = u − (that multiple of K). An odd
# multiple of K is a quarter period away from a zero of sn: with sn(K + v) = cd v,
# cn(K + v) = −k′ sd v, dn(K + v) = k′ nd v, the reduced argument is v ± K, and cn vanishes
# exactly at the turning point.
function _phase_jacobi(q, scale, K, J, kp, polar)
    j = round(q / pi)
    v = scale * (q - j * pi)
    nK = polar ? 2j + 1 : j
    sv, cv, dv = _jacobi_sncndn(v, J)
    if isodd(Int(rem(nK, 2)))
        ε = v <= 0 ? one(v) : -one(v)
        return (ε * cv / dv, -ε * kp * sv / dv, kp / dv), v + ε * K, isodd(Int(rem((nK - ε) / 2, 2)))
    end
    return (sv, cv, dv), v, isodd(Int(rem(nK / 2, 2)))
end

# The radial motion between the turning points r2 (q = 0) and r1 (q = π).
struct _RadialPhaseModel{T}
    roots::NTuple{4,T}
    k::T; m1::T; K::T; scale::T; kp::T
    J::_JacobiParameter{T}
    h::T; n1::T; P::T; D::T            # characteristic hr, 1 − hr, P(hr|k), D(k)
    # the horizon terms as Σ wᵢ F(ρᵢ) over the nodes of the divided difference between r₊ and r₋
    nh::Vector{T}; n1h::Vector{T}; Ph::Vector{T}; wt::Vector{T}; wϕ::Vector{T}
    tfactor::T; ϕfactor::T             # −E/√((1 − E²)(r1 − r3)(r2 − r4)), 2aE(r2 − r3)/(…)
    ϒ::T
end

function _radial_phase_model(a, E, Lz, ν, roots, ϒr)
    r1, r2, r3, r4 = roots
    T = typeof(E)
    k = ((r1 - r2) / (r1 - r3)) * ((r3 - r4) / (r2 - r4))
    m1 = (r2 - r3) * (r1 - r4) / ((r1 - r3) * (r2 - r4))
    K = _ellip_k(m1)
    J = k <= 1 ? _JacobiParameter{T}(k, one(T), one(T), _landen(k, m1)) : _jacobi_parameter(k)
    h = (r1 - r2) / (r1 - r3)
    n1 = (r2 - r3) / (r1 - r3)
    pole(n1) = _carlson_rj(zero(T), m1, one(T), n1) / 3
    rp = 1 + sqrt(1 - a^2)
    nodes, weights = _horizon_divdiff_rule(rp, a^2 / rp)
    nh = [(r1 - r2) * (r3 - ρ) / ((r1 - r3) * (r2 - ρ)) for ρ in nodes]
    n1h = [(r2 - r3) * (r1 - ρ) / ((r1 - r3) * (r2 - ρ)) for ρ in nodes]
    aLE = a * Lz / E
    wt = [w * (-2a^2 + ρ * (4 - aLE)) / ((r2 - ρ) * (r3 - ρ)) for (ρ, w) in zip(nodes, weights)]
    wϕ = [w * (2ρ - aLE) / ((r2 - ρ) * (r3 - ρ)) for (ρ, w) in zip(nodes, weights)]
    sqrtfac = sqrt(max(zero(T), ν * (r1 - r3) * (r2 - r4)))
    return _RadialPhaseModel{T}(roots, k, m1, K, K / pi, sqrt(m1), J, h, n1, pole(n1),
        _ellip_d_complete(m1), nh, n1h, pole.(n1h), wt, wϕ, -E / sqrtfac,
        2a * E * (r2 - r3) / sqrtfac, ϒr)
end

# The polar motion z = z₋ sn(K(2q/π + 1)), z = z₋ at q = 0.
struct _PolarPhaseModel{T}
    zp::T; zm::T; x::T
    k::T; m1::T; K::T; scale::T; kp::T
    J::_JacobiParameter{T}
    n::T; n1::T; P::T; D::T            # characteristic 1 − x², x², P(1 − x²|k) (L_z ≠ 0), D(k)
    tfactor::T; ϕfactor::T             # E a² z₋²/z₊, −L_z (1 − x²)/z₊
    side::T                            # the side of L_z = 0, as in ϒϕ
    axis::Bool                         # L_z = 0
    ϒ::T
end

function _polar_phase_model(a, x, E, Lz, ν, L, ϒθ)
    T = typeof(E)
    zp, zm = sqrt(a^2 * ν + L^2), sqrt(1 - x^2)
    k = a^2 * ν * (zm / zp)^2
    m1 = 1 - k
    K = _ellip_k(m1)
    n, n1 = 1 - x^2, x^2
    axis = iszero(Lz)
    P = axis ? zero(T) : _carlson_rj(zero(T), m1, one(T), n1) / 3
    return _PolarPhaseModel{T}(zp, zm, x, k, m1, K, 2K / pi, sqrt(m1), _jacobi_parameter(k),
        n, n1, P, _ellip_d_complete(m1), E * a^2 * zm^2 / zp, -Lz * n / zp, _polar_side(x),
        axis, ϒθ)
end

"""
    KerrGeoOrbitPhases{T}

The orbit of [`kerr_geo_orbit`](@ref) as a record in the floating-point type `T`, evaluated by
its radial phase ``q_r`` with [`kerr_geo_radial_motion`](@ref) and by its polar phase
``q_\\theta`` with [`kerr_geo_polar_motion`](@ref). Built by [`kerr_geo_orbit_phases`](@ref).

Fields: `a`, `p`, `e`, `x`; the constants of motion `E`, `Lz`, `Q`; the Mino frequencies
`Upsilon_t`, `Upsilon_r`, `Upsilon_theta`, `Upsilon_phi` (those of
[`kerr_geo_frequencies`](@ref), with `Upsilon_r = 0` for a constant-radius Critical orbit);
`radial_roots` `(r1, r2, r3, r4)` and `polar_roots` `(zp, zm)`; `precision`, the `BigFloat`
precision the record evaluates at (0 for other types). The remaining fields are the
precomputed elliptic data of the evaluation.

The record is immutable: one record can be evaluated from several tasks at once.
"""
struct KerrGeoOrbitPhases{T}
    a::T; p::T; e::T; x::T
    E::T; Lz::T; Q::T
    Upsilon_t::T; Upsilon_r::T; Upsilon_theta::T; Upsilon_phi::T
    radial_roots::NTuple{4,T}
    polar_roots::NTuple{2,T}
    precision::Int
    radial::_RadialPhaseModel{T}
    polar::_PolarPhaseModel{T}
end

"""
    kerr_geo_orbit_phases(a, p, e, x; precision=nothing) -> KerrGeoOrbitPhases

The orbit that [`kerr_geo_orbit`](@ref) builds for `(a, p, e, x)`, a stable orbit or a
constant-radius Critical orbit with E < 1, as a [`KerrGeoOrbitPhases`](@ref) record evaluated
by its phases. The phases are those of `kerr_geo_orbit` without initial phases: ``q_r = 0`` at
the periapsis ``r_2`` and ``q_\\theta = 0`` at the polar turning point ``z = z_-``, and
``q = \\Upsilon \\lambda`` along the orbit. Any other orbit throws the `ArgumentError` of
`kerr_geo_orbit`, and [`kerr_geo_orbit_supported`](@ref) tells which orbits these are.

The record is computed in the floating-point type of `(a, p, e, x)`; `precision = p` converts
them to `BigFloat` of `p` bits, and the record then evaluates at that precision.
"""
function kerr_geo_orbit_phases(a::Real, p::Real, e::Real, x::Real; precision=nothing)
    precision === nothing || return setprecision(BigFloat, precision) do
        kerr_geo_orbit_phases(BigFloat(a), BigFloat(p), BigFloat(e), BigFloat(x))
    end
    T = _float_type(a, p, e, x)
    a, p, e, x = T(a), T(p), T(e), T(x)
    e > 1 && throw(ArgumentError("kerr_geo_orbit: e = $e > 1 is an unbound orbit; the bound-orbit " *
        "formulas hold for e ≤ 1. Use kerr_geodesic(a, p, e, x).Scatter."))
    prec = _input_precision(a, p, e, x)
    return _with_precision(T, prec) do
        meta = _apex_orbit_gate(a, p, e, x)
        F = _kerr_geo_frequencies(a, p, e, x, "Mino")
        ϒr = meta.family == "Critical" ? zero(T) : T(F["ϒr"])
        c = _apex_constants(a, p, e, x)
        roots = kerr_geo_radial_roots(a, p, e, x; En=c.E, Lz=c.Lz, Q=c.Q)
        radial = _radial_phase_model(a, c.E, c.Lz, c.ν, roots, ϒr)
        polar = _polar_phase_model(a, x, c.E, c.Lz, c.ν, c.L, T(F["ϒθ"]))
        KerrGeoOrbitPhases{T}(a, p, e, x, c.E, c.Lz, c.Q, T(F["ϒt"]), ϒr, T(F["ϒθ"]),
            T(F["ϒϕ"]), roots, (polar.zp, polar.zm), prec, radial, polar)
    end
end

Base.show(io::IO, o::KerrGeoOrbitPhases{T}) where {T} =
    print(io, "KerrGeoOrbitPhases{", T, "}(a=", o.a, ", p=", o.p, ", e=", o.e, ", x=", o.x, ")")

function _radial_motion(M::_RadialPhaseModel, q)
    r1, r2, r3, r4 = M.roots
    (s, c, d), y, _ = _phase_jacobi(q, M.scale, M.K, M.J, M.kp, false)
    s2, c2, d2 = s^2, c^2, d^2
    yK = y / M.K
    # r − r2 = (r1 − r2)(r2 − r3) sn²/((r2 − r3) + (r1 − r2) cn²)
    den = (r2 - r3) + (r1 - r2) * c2
    r = r2 + (r1 - r2) * (r2 - r3) * s2 / den
    drdq = 2 * (r1 - r2) * (r1 - r3) * (r2 - r3) * M.scale * s * c * d / den^2
    # Π(h|k) y/K − Π(h; am y|k) and its y derivative, with 1 − h sn² = cn² + (1 − h) sn²
    I(h, n1, P) = h * (P * yK - _ellip_pole(s, c, M.m1, n1))
    dI(h, n1, P) = h * (P / M.K - s2 / (c2 + n1 * s2))
    Ir, dIr = I(M.h, M.n1, M.P), dI(M.h, M.n1, M.P)
    Ht = Hϕ = dHt = dHϕ = zero(r)
    for i in eachindex(M.nh)
        Ii, dIi = I(M.nh[i], M.n1h[i], M.Ph[i]), dI(M.nh[i], M.n1h[i], M.Ph[i])
        Ht += M.wt[i] * Ii; Hϕ += M.wϕ[i] * Ii
        dHt += M.wt[i] * dIi; dHϕ += M.wϕ[i] * dIi
    end
    g = c2 + M.n1 * s2
    A = (r2 - r3) * (4 + r1 + r2 + r3 + r4)
    B = (r1 - r3) * (r2 - r4)
    term = B * (M.h * s * c * d / g - M.k * (M.D * yK - _ellip_d(s, c, M.m1)))
    dterm = B * (M.h * (c2 * d2 - s2 * d2 - M.k * s2 * c2) / g + 2 * M.h^2 * s2 * c2 * d2 / g^2 -
        M.k * (M.D / M.K - s2))
    return (r = r, dr_dlambda = M.ϒ * drdq,
        delta_t = M.tfactor * (A * Ir - 4 * (r2 - r3) * Ht + term),
        delta_phi = M.ϕfactor * Hϕ,
        ddelta_t_dq = M.tfactor * M.scale * (A * dIr - 4 * (r2 - r3) * dHt + dterm),
        ddelta_phi_dq = M.ϕfactor * M.scale * dHϕ)
end

function _polar_motion(M::_PolarPhaseModel, q)
    (s, c, d), y, odd = _phase_jacobi(q, M.scale, M.K, M.J, M.kp, true)
    σ = odd ? -one(s) : one(s)
    s2 = s^2
    yK = y / M.K
    z = M.zm * σ * s
    dzdq = M.zm * σ * c * d * M.scale
    # sin²θ = 1 − z₋² sn² = cn² + x² sn²; at a turning point dθ/dλ = 0, the mean of its two
    # sides where the orbit passes over the axis (x = 0)
    dθdq = iszero(c) ? zero(s) : -M.zm * d * M.scale * (σ * c) / hypot(c, M.x * s)
    if M.axis
        # L_z = 0: φ steps by side·π at every pass over the axis, q = jπ, and takes the mean
        # of its two sides on the pass, as the polar part of φ in kerr_geo_orbit
        δ = q - round(q / pi) * pi
        Δϕ, dΔϕ = M.side * (sign(δ) * pi / 2 - δ), -M.side
    else
        Δϕ = M.ϕfactor * (M.P * yK - _ellip_pole(s, c, M.m1, M.n1))
        dΔϕ = M.ϕfactor * M.scale * (M.P / M.K - s2 / (c^2 + M.n1 * s2))
    end
    return (z = z, theta = acos(z), dz_dlambda = M.ϒ * dzdq, dtheta_dlambda = M.ϒ * dθdq,
        delta_t = M.tfactor * (_ellip_d(s, c, M.m1) - M.D * yK),
        delta_phi = Δϕ,
        ddelta_t_dq = M.tfactor * M.scale * (s2 - M.D / M.K),
        ddelta_phi_dq = dΔϕ)
end

const _RADIAL_MOTION_FIELDS = (:r, :dr_dlambda, :delta_t, :delta_phi, :ddelta_t_dq, :ddelta_phi_dq)
const _POLAR_MOTION_FIELDS = (:z, :theta, :dz_dlambda, :dtheta_dlambda, :delta_t, :delta_phi,
    :ddelta_t_dq, :ddelta_phi_dq)

"""
    kerr_geo_radial_motion(orbit::KerrGeoOrbitPhases, qr)
    kerr_geo_radial_motion(orbit, qr::AbstractVector)
    kerr_geo_radial_motion!(out, orbit, qr::AbstractVector)

The radial motion of `orbit` at the radial phase `qr` (``q_r = 0`` at the periapsis ``r_2``,
``q_r = \\pi`` at the apoapsis ``r_1``), as the `NamedTuple`

- `r`, `dr_dlambda` = ``dr/d\\lambda = \\pm\\sqrt{R(r)}``;
- `delta_t`, `delta_phi`: the radial parts ``\\Delta t_r``, ``\\Delta\\phi_r`` of t and φ,
  zero at ``q_r = 0`` (`"CrossFunction"[1]` and `[3]` of [`kerr_geo_orbit`](@ref));
- `ddelta_t_dq`, `ddelta_phi_dq`: their ``q_r`` derivatives (`"DerivativesCrossFunction"`).

Along the orbit ``t = \\Upsilon_t \\lambda + \\Delta t_r(q_r) + \\Delta t_\\theta(q_\\theta)`` and
likewise φ, with ``q = \\Upsilon \\lambda``. At ``e = 0`` the radius is constant and the radial
velocity and parts vanish identically. For a vector of phases the method returns a
`NamedTuple` of vectors, and `kerr_geo_radial_motion!` fills a preallocated one with the same
field names. In `Float64` the scalar call does not allocate.
"""
kerr_geo_radial_motion(o::KerrGeoOrbitPhases{T}, q::Real) where {T} =
    _with_precision(() -> _radial_motion(o.radial, T(q)), T, o.precision)

"""
    kerr_geo_polar_motion(orbit::KerrGeoOrbitPhases, qθ)
    kerr_geo_polar_motion(orbit, qθ::AbstractVector)
    kerr_geo_polar_motion!(out, orbit, qθ::AbstractVector)

The polar motion of `orbit` at the polar phase `qθ` (``q_\\theta = 0`` at the turning point
``z = z_-``, ``q_\\theta = \\pi`` at ``z = -z_-``), as the `NamedTuple`

- `z` = ``\\cos\\theta``, `theta`, `dz_dlambda` = ``\\pm\\sqrt{\\Theta(z)}`` with
  ``\\Theta(z) = Q - z^2(a^2(1 - E^2)(1 - z^2) + L_z^2 + Q)``, and `dtheta_dlambda`;
- `delta_t`, `delta_phi`: the polar parts ``\\Delta t_\\theta``, ``\\Delta\\phi_\\theta`` of t and
  φ, zero at ``q_\\theta = 0`` (`"CrossFunction"[2]` and `[4]` of [`kerr_geo_orbit`](@ref));
- `ddelta_t_dq`, `ddelta_phi_dq`: their ``q_\\theta`` derivatives.

At ``|x| = 1`` the orbit is equatorial, θ = π/2, and the polar velocity and parts vanish
identically. At `x = 0` (``L_z = 0``) every value is finite: the orbit passes over the axis at
``q_\\theta = j\\pi``, where `dtheta_dlambda` changes sign and `delta_phi` steps by π (by −π
for `x = -0.0`), the limit ``x \\to 0^\\pm`` of `kerr_geo_orbit`. On the pass itself both take
the mean of their two sides: `dtheta_dlambda` = 0, and `delta_phi` lies halfway, so that the
side ``q_\\theta \\to j\\pi^\\pm`` is `delta_phi ± π/2` (`delta_phi ∓ π/2` for `x = -0.0`).
`ddelta_phi_dq` is ``\\mp 1`` there, without the step. Vectors of phases as in
[`kerr_geo_radial_motion`](@ref).
"""
kerr_geo_polar_motion(o::KerrGeoOrbitPhases{T}, q::Real) where {T} =
    _with_precision(() -> _polar_motion(o.polar, T(q)), T, o.precision)

function _motion!(f, names, out::NamedTuple, model, q::AbstractVector)
    keys(out) == names || throw(ArgumentError("the output needs the fields $names"))
    all(v -> axes(v) == axes(q), values(out)) ||
        throw(DimensionMismatch("each output vector needs the axes of the phases"))
    for i in eachindex(q)
        v = f(model, q[i])
        map((dst, val) -> (dst[i] = val), values(out), values(v))
    end
    return out
end

"""
    kerr_geo_radial_motion!(out, orbit::KerrGeoOrbitPhases, qr::AbstractVector)

Fill the `NamedTuple` of vectors `out`, with the fields of [`kerr_geo_radial_motion`](@ref)
in its order and the axes of `qr`, with the radial motion at the phases `qr`. Returns `out`;
in `Float64` it does not allocate.
"""
kerr_geo_radial_motion!(out::NamedTuple, o::KerrGeoOrbitPhases{T}, q::AbstractVector) where {T} =
    _with_precision(() -> _motion!((M, qi) -> _radial_motion(M, T(qi)), _RADIAL_MOTION_FIELDS,
        out, o.radial, q), T, o.precision)
"""
    kerr_geo_polar_motion!(out, orbit::KerrGeoOrbitPhases, qθ::AbstractVector)

Fill the `NamedTuple` of vectors `out`, with the fields of [`kerr_geo_polar_motion`](@ref) in
its order and the axes of `qθ`, with the polar motion at the phases `qθ`. Returns `out`; in
`Float64` it does not allocate.
"""
kerr_geo_polar_motion!(out::NamedTuple, o::KerrGeoOrbitPhases{T}, q::AbstractVector) where {T} =
    _with_precision(() -> _motion!((M, qi) -> _polar_motion(M, T(qi)), _POLAR_MOTION_FIELDS,
        out, o.polar, q), T, o.precision)

_motion_arrays(::Type{T}, names, q) where {T} =
    NamedTuple{names}(ntuple(_ -> similar(q, T), length(names)))
kerr_geo_radial_motion(o::KerrGeoOrbitPhases{T}, q::AbstractVector) where {T} =
    kerr_geo_radial_motion!(_motion_arrays(T, _RADIAL_MOTION_FIELDS, q), o, q)
kerr_geo_polar_motion(o::KerrGeoOrbitPhases{T}, q::AbstractVector) where {T} =
    kerr_geo_polar_motion!(_motion_arrays(T, _POLAR_MOTION_FIELDS, q), o, q)
