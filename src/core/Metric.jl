# Kerr metric functions shared by everything: horizons, energy regime, radial and polar
# potentials and their roots, and the tortoise coordinate. The degeneracy tolerances it uses
# are defined in core/Degeneracy.jl.

"""
    kerr_delta(a, r)

``\\Delta(r) = r^2 - 2r + a^2 = (r - r_+)(r - r_-)``, evaluated in the product form, which keeps its
relative precision next to a horizon.
"""
function kerr_delta(a::Real, r::Real)
    # r^2 - 2r + a^2 cancels to an O(eps) absolute error while Δ -> 0 at a horizon
    s = sqrt(max((1 - a) * (1 + a), zero(float(a))))
    return (r - 1 - s) * (r - 1 + s)
end

_rplus(a) = 1 + sqrt(1 - a^2)
_rminus(a) = 1 - sqrt(1 - a^2)

"""
    kerr_horizons(a)

Horizon radii ``r_\\pm = 1 \\pm \\sqrt{1 - a^2}`` for ``|a| \\leq 1``
(``G = c = M = 1``), returned as `(rplus, rminus)`.
"""
function kerr_horizons(a::Real; atol::Real=_classification_atol(float(typeof(a))))
    abs(a) <= 1 + atol || throw(DomainError(a, "Kerr spin must satisfy |a|<=1."))
    spin = clamp(float(a), -1, 1)
    return (rplus=_rplus(spin), rminus=_rminus(spin))
end

"""
    kerr_metric_limit(a; atol=64eps(), near_extremal_threshold=1e-6)

`:schwarzschild` (|a| ≤ `atol`), `:extremal` (||a| − 1| ≤ `atol`), `:near_extremal`
(1 − |a| ≤ `near_extremal_threshold`) or `:subextremal`, tested in that order.
"""
function kerr_metric_limit(a::Real;
        atol::Real=_classification_atol(float(typeof(a))),
        near_extremal_threshold::Real=1.0e-6)
    abs(a) <= 1 + atol || throw(DomainError(a, "Kerr spin must satisfy |a|<=1."))
    abs(a) <= atol && return :schwarzschild
    abs(abs(a) - 1) <= atol && return :extremal
    1 - abs(a) <= near_extremal_threshold && return :near_extremal
    return :subextremal
end

"""
    _e2m1(E)

E² − 1 as (E − 1)(E + 1), accurate to a few ulps for every E (E − 1 is exact near 1). The
quartic coefficient of R and the polar constant a²(1 − E²) use it: E² formed first rounds
away the δ² of E = 1 + δ, a relative error of up to ~5e-9 (at |δ| ≈ 1e-8).
"""
_e2m1(energy) = (energy - 1) * (energy + 1)

"""
    kerr_energy_regime(E; atol=0, rtol=0)

`:elliptic`, `:parabolic` or `:hyperbolic` by the sign of ``E^2 - 1``, the ``r^4`` coefficient of ``R``
(negative, zero or positive), for either sign of E; |E| within `atol + rtol·max(1, |E|)` of 1
is `:parabolic`. The sign of E itself is `kerr_energy_sign`.
"""
function kerr_energy_regime(energy::Real;
        atol::Real=DEFAULT_ENERGY_ATOL,
        rtol::Real=DEFAULT_ENERGY_RTOL)
    tolerance = atol + rtol * max(1.0, abs(float(energy)))
    distance = abs(energy) - 1                  # the sign of E² − 1
    distance < -tolerance && return :elliptic
    distance > tolerance && return :hyperbolic
    return :parabolic
end

"""
    kerr_energy_sign(E)

`+1` for ``E \\geq 0`` and `-1` for ``E < 0``. Future-directed motion with ``E < 0`` exists only inside the
ergoregion: the Trapped class.
"""
kerr_energy_sign(energy::Real) = energy < 0 ? -1 : 1

"""
    kerr_radial_momentum(a, E, Lz, r)

``P(r) = E(r^2 + a^2) - aL_z``. It enters the radial potential,
``R = P^2 - \\Delta[r^2 + (L_z - aE)^2 + Q]``, and the Mino-time rates
``dt/d\\lambda`` and ``d\\phi/d\\lambda``; on the horizon ``R(r_+) = P(r_+)^2``.
"""
kerr_radial_momentum(a::Real, energy::Real, lz::Real, r::Real) =
    energy * (r - 1) * (r + 1) + (energy * (1 + a^2) - a * lz)   # r² − 1 accurate near r = 1

"""
    kerr_axis_carter_q(a, E)

Carter constant ``Q = a^2(1 - E^2)`` of a timelike geodesic along the spin axis (``L_z = 0``).
"""
kerr_axis_carter_q(a::Real, energy::Real) = -a^2 * _e2m1(energy)

"""
    kerr_axis_radial_potential(a, E, r)

Radial potential on the spin axis (``L_z = 0``, ``Q = a^2(1 - E^2)``) in factored form,
``R = (r^2 + a^2)[E^2(r^2 + a^2) - \\Delta]``.
"""
function kerr_axis_radial_potential(a::Real, energy::Real, r::Real)
    sigma = r^2 + a^2
    return sigma * (energy^2 * sigma - kerr_delta(a, r))
end

"""
    kerr_radial_coefficients(a, E, Lz, Q; energy_atol=0, energy_rtol=0)

Coefficients ``(c_0, c_1, c_2, c_3, c_4)`` of ``R(r) = \\sum_{k=0}^4 c_k r^k``,
in ascending order. ``c_4 = E^2 - 1`` is
set to zero when `kerr_energy_regime(E; atol=energy_atol, rtol=energy_rtol)` is
`:parabolic`, so ``R`` is built as a cubic; with the default zero tolerances this is ``E^2 = 1``
exactly.
"""
function kerr_radial_coefficients(a::Real, energy::Real, lz::Real, q::Real;
        energy_atol::Real=DEFAULT_ENERGY_ATOL,
        energy_rtol::Real=DEFAULT_ENERGY_RTOL)
    T = _float_type(a, energy, lz, q)
    regime = kerr_energy_regime(energy; atol=energy_atol, rtol=energy_rtol)
    c4 = regime === :parabolic ? zero(T) : T(_e2m1(energy))
    return (
        T(-a^2 * q),
        T(2 * (a * energy - lz)^2 + 2q),
        T(-(q + lz^2 - a^2 * _e2m1(energy))),
        T(2),
        c4,
    )
end

# R(1+x), formed from P(1) and Delta(1). Unlike translating the expanded
# r-polynomial, this preserves P(1)^2 when two roots approach the extremal horizon.
function _radial_shifted_coefficients(a,energy,lz,q)
    p=kerr_radial_momentum(a,energy,lz,one(_float_type(a,energy,lz,q)))
    d=(a-1)*(a+1)
    k=1+(lz-a*energy)^2+q
    return (p^2-d*k,4energy*p-2d,4energy^2+2energy*p-k-d,
        4energy^2-2,_e2m1(energy))
end

"""
    kerr_radial_polynomial(a, E, Lz, Q; energy_atol=0, energy_rtol=0)

``R(r)`` as a `Polynomial` (Polynomials.jl) built from `kerr_radial_coefficients`; it is a
cubic when ``E^2 = 1``.
"""
function kerr_radial_polynomial(a::Real, energy::Real, lz::Real, q::Real; kwargs...)
    coefficients = collect(kerr_radial_coefficients(a, energy, lz, q; kwargs...))
    while length(coefficients) > 1 && iszero(coefficients[end])
        pop!(coefficients)
    end
    return Polynomial(coefficients)
end

"""
    kerr_radial_potential(a, E, Lz, Q, r)

Radial potential of a timelike Kerr geodesic with Carter constant ``Q``,
``R = [E(r^2 + a^2) - aL_z]^2 - \\Delta[r^2 + (L_z - aE)^2 + Q]``. It is summed from its coefficients
(`kerr_radial_coefficients`), which keeps full relative precision at large r, where the
factored form would cancel two terms of size ``E^2r^4``.
"""
function kerr_radial_potential(a::Real, energy::Real, lz::Real, q::Real, r::Real)
    c0, c1, c2, c3, c4 = kerr_radial_coefficients(a, energy, lz, q)
    return c0 + r * (c1 + r * (c2 + r * (c3 + r * c4)))
end

"""
    kerr_radial_derivatives(a, E, Lz, Q, r; energy_atol=0, energy_rtol=0)

R and its first four derivatives at `r`, as the NamedTuple `(R, R1, R2, R3, R4)`.
"""
function kerr_radial_derivatives(a::Real, energy::Real, lz::Real, q::Real, r::Real;
        kwargs...)
    c0, c1, c2, c3, c4 = kerr_radial_coefficients(a, energy, lz, q; kwargs...)
    return (
        R=c0 + c1 * r + c2 * r^2 + c3 * r^3 + c4 * r^4,
        R1=c1 + 2 * c2 * r + 3 * c3 * r^2 + 4 * c4 * r^3,
        R2=2 * c2 + 6 * c3 * r + 12 * c4 * r^2,
        R3=6 * c3 + 24 * c4 * r,
        R4=24 * c4,
    )
end

"""
    _polish_root(coefficients, z; order=0)

Newton's method on the radial polynomial with the given coefficients (constant term first), or
on its `order`-th derivative (order 1 for a double root), from the estimate `z` (real or
complex): at most eight steps, stopping when the step is within 4 eps of `z` or when the
derivative is below rounding relative to the size of its terms.
"""
function _polish_root(coefficients, z; order::Int=0, evaluator=nothing)
    T = real(float(typeof(z)))
    c = coefficients
    for _ in 1:order
        c = ntuple(i -> i * c[i + 1], length(c) - 1)
    end
    dc = ntuple(i -> i * c[i + 1], length(c) - 1)
    for _ in 1:8
        az = abs(z)
        scale = sum(abs(dc[i]) * az^(i - 1) for i in eachindex(dc))
        value,dp = evaluator===nothing ? (evalpoly(z,c),evalpoly(z,dc)) : evaluator(z)
        abs(dp) <= eps(T) * scale && break
        step = value / dp
        z -= step
        abs(step) <= 4 * eps(T) * max(1.0, abs(z)) && break
    end
    return z
end

# Error-free transforms keep coefficient formation and Newton residuals accurate
# when simple roots nearly coincide. The Newton iteration itself is unchanged.
function _two_sum(a,b)
    s=a+b
    v=s-a
    return s,(a-(s-v))+(b-v)
end
_two_product(a::Real,b::Real) = (a*b,fma(a,b,-a*b))
function _two_product(a::Complex,b::Complex)
    ac,eac=_two_product(real(a),real(b)); bd,ebd=_two_product(imag(a),imag(b))
    ad,ead=_two_product(real(a),imag(b)); bc,ebc=_two_product(imag(a),real(b))
    re,er=_two_sum(ac,-bd); im,ei=_two_sum(ad,bc)
    return complex(re,im),complex(er+eac-ebd,ei+ead+ebc)
end
_two_product(a::Real,b::Complex) = _two_product(complex(a),b)
_two_product(a::Complex,b::Real) = _two_product(a,complex(b))
_wide(x)=(x,zero(x))
function _wide_add(a,b)
    s,e=_two_sum(a[1],b[1])
    return _two_sum(s,e+a[2]+b[2])
end
_wide_neg(a)=(-a[1],-a[2])
function _wide_sqrt(x)
    root = sqrt(x[1])
    return _two_sum(root, (fma(-root, root, x[1]) + x[2]) / (2root))
end
function _wide_div(a,b)
    quotient = a[1] / b[1]
    rest = _wide_sub(a, _wide_mul(_wide(quotient), b))
    return _two_sum(quotient, (rest[1] + rest[2]) / b[1])
end
_wide_sub(a,b)=_wide_add(a,_wide_neg(b))
function _wide_mul(a,b)
    p,e=_two_product(a[1],b[1])
    return _two_sum(p,e+a[1]*b[2]+a[2]*b[1]+a[2]*b[2])
end
function _wide_evalpoly(x,c)
    z=_wide(x); v=c[end]
    for i in length(c)-1:-1:1
        v=_wide_add(_wide_mul(v,z),c[i])
    end
    return v[1]+v[2]
end

# The coefficients of R in double-double, from the exact inputs.
function _wide_radial_coefficients(a,E,L,Q)
    aa,ee,ll,qq=_wide.(float.((a,E,L,Q)))
    a2=_wide_mul(aa,aa)
    lead=_wide_mul(_wide_sub(ee,_wide(1.0)),_wide_add(ee,_wide(1.0)))
    u=_wide_sub(ll,_wide_mul(aa,ee))
    return (_wide_neg(_wide_mul(a2,qq)),
       _wide_mul(_wide(2.0),_wide_add(_wide_mul(u,u),qq)),
       _wide_sub(_wide_mul(a2,lead),_wide_add(_wide_mul(ll,ll),qq)),
       _wide(2.0),lead)
end

# coefficients of R^(k): j!/(j − k)! c_j
_wide_derivative_coefficients(c,k) =
    ntuple(i->_wide_mul(_wide(float(factorial(i+k-1)÷factorial(i-1))),c[i+k]),length(c)-k)

function _radial_root_evaluator(a,E,L,Q)
    c=_wide_radial_coefficients(a,E,L,Q)
    dc=_wide_derivative_coefficients(c,1)
    return r->(_wide_evalpoly(r,c),_wide_evalpoly(r,dc))
end

# r₊ = 1 + √(1 − a²) as a double-double pair, and P(r₊) = 2E r₊ − aLz (r₊² + a² = 2r₊)
# from it: P(r₊) cancels when the horizon is nearly a root, and the rounding of r₊ alone would
# leave ~eps·E of it.
function _wide_rplus(a)
    x = _wide_mul(_two_sum(1.0, -float(a)), _two_sum(1.0, float(a)))      # 1 − a², exactly paired
    root = sqrt(x[1])
    root_low = iszero(root) ? zero(root) : (fma(-root, root, x[1]) + x[2]) / (2root)   # |a| = 1: r₊ = 1
    return _wide_add(_wide(1.0), _two_sum(root, root_low)), 2 * (root + root_low)
end

function _wide_horizon_momentum(a, E, L)
    rplus, _ = _wide_rplus(a)
    p = _wide_sub(_wide_mul(_wide(2.0 * E), rplus), _wide_mul(_wide(float(a)), _wide(float(L))))
    return p[1] + p[2]
end

"""
    _horizon_gap(a, E, Lz, Q, r)

r − r₊ for a simple root r of R, to relative rounding even when the gap is far below ulp(r₊)
(a turning point next to a nearly-root horizon, P(r₊) → 0). Newton, from the rounded root on
either side of r₊, on R about the horizon in double-double,
    R(r₊ + δ) = P₊² + (4E r₊P₊ − dK₊) δ + (2EP₊ + 4E²r₊² − K₊ − 2d r₊) δ² + (4E²r₊ − 2r₊ − d) δ³
                + (E² − 1) δ⁴,   d = r₊ − r₋,
which uses Δ(r₊) = 0 and r₊² + a² = 2r₊ exactly; the coefficients are formed from the exact
inputs, so a root far from r₊ is as accurate as the double-double root refinement.
"""
function _horizon_gap(a, E, L, Q, r)
    h, _ = _wide_rplus(a)
    d = _wide_mul(_wide(2.0), _wide_sub(h, _wide(1.0)))          # r₊ − r₋ = 2(r₊ − 1)
    ee = _wide(float(E))
    p = _wide_sub(_wide_mul(_wide(2.0), _wide_mul(ee, h)), _wide_mul(_wide(float(a)), _wide(float(L))))
    u = _wide_sub(_wide(float(L)), _wide_mul(_wide(float(a)), ee))
    k = _wide_add(_wide_add(_wide_mul(h, h), _wide_mul(u, u)), _wide(float(Q)))
    times(n, x) = _wide_mul(_wide(float(n)), x)
    eh = _wide_mul(ee, h); e2 = _wide_mul(ee, ee)
    c = (_wide_mul(p, p),
         _wide_sub(times(4, _wide_mul(eh, p)), _wide_mul(d, k)),
         _wide_sub(_wide_add(times(4, _wide_mul(eh, eh)), times(2, _wide_mul(ee, p))),
                   _wide_add(k, times(2, _wide_mul(d, h)))),
         _wide_sub(times(4, _wide_mul(eh, ee)), _wide_add(times(2, h), d)),
         _wide_mul(_wide_sub(ee, _wide(1.0)), _wide_add(ee, _wide(1.0))))
    dc = (c[2], times(2, c[3]), times(3, c[4]), times(4, c[5]))
    delta = r - (h[1] + h[2])
    for _ in 1:64
        step = _wide_evalpoly(delta, c) / _wide_evalpoly(delta, dc)
        delta -= step
        abs(step) <= 4eps(delta) && break
    end
    return delta
end

# All roots of R refined together (Aberth–Ehrlich) with the double-double residual of
# `evaluator`. The correction N/(1 − N Σ_{j≠i} 1/(z_i − z_j)), N = R/R′, keeps the estimates
# apart, so roots that cluster (near r = 0 for small spin and small Q + (Lz − aE)², near r = 1
# close to extremality) converge to the roots of the exact-coefficient polynomial; Newton on
# one root at a time slides into the cluster instead. Started from the companion roots of
# R(1 + x); a double root is approached by two estimates, which the repeated-root
# classification then takes together.
function _refine_roots(evaluator, z::Vector{Complex{T}}) where {T}
    n = length(z)
    # Nearly equal roots can need more than 64 sweeps to separate.
    for _ in 1:1000
        largest = zero(T)
        for i in 1:n
            value, slope = evaluator(z[i])
            iszero(value) && continue
            newton = value / slope
            pull = sum((inv(z[i] - z[j]) for j in 1:n if j != i); init=zero(Complex{T}))
            step = newton / (1 - newton * pull)
            isfinite(step) || continue
            z[i] -= step
            largest = max(largest, abs(step) / max(abs(z[i]), floatmin(T)))
        end
        largest <= 4eps(T) && break
    end
    return z
end

# Roots of a real polynomial (ascending coefficients) in the floating-point type T of the
# coefficients: the companion-matrix roots of the coefficients rounded to Float64, each
# polished by Newton's method in T. A step is taken only while it exceeds the rounding of the
# root, so a root that is already accurate to eps(T) is returned as found.
function _polynomial_roots(coefficients)
    T = float(eltype(coefficients))
    c = collect(T, coefficients)
    while length(c) > 1 && iszero(c[end])
        pop!(c)
    end
    dc = [k * c[k + 1] for k in 1:length(c) - 1]
    estimates = roots(Polynomial(Float64.(c)))
    return map(estimates) do z0
        z = Complex{T}(z0)
        for _ in 1:8 + 2 * ceil(Int, log2(precision(T) / 53))
            step = evalpoly(z, c) / evalpoly(z, dc)
            (isfinite(step) && abs(step) > 4eps(T) * abs(z)) || break
            z -= step
        end
        z
    end
end

function _derivative_scales(coefficients, r)
    c0, c1, c2, c3, c4 = coefficients
    ar = abs(float(r))
    return (
        abs(c0) + abs(c1)*ar + abs(c2)*ar^2 + abs(c3)*ar^3 + abs(c4)*ar^4,
        abs(c1) + 2 * abs(c2) * ar + 3 * abs(c3) * ar^2 + 4 * abs(c4) * ar^3,
        2 * abs(c2) + 6 * abs(c3) * ar + 12 * abs(c4) * ar^2,
        6 * abs(c3) + 24 * abs(c4) * ar,
        24 * abs(c4),
    )
end

"""
    kerr_root_multiplicity_at(a, E, Lz, Q, r; atol=1e-12, rtol=1e-12,
                              energy_atol=0, energy_rtol=0)

Multiplicity of the radius `r` as a root of R: the number of consecutive values R, R′, R″,
R‴, R⁗ at `r`, starting from R, that vanish within `atol + rtol·max(1, s)`, where `s` is the
sum of the magnitudes of that derivative's terms (0 when R(r) ≠ 0). The default tolerances
are 1e-12 in Float64, carried to other floating-point types by `_tol`.
"""
function kerr_root_multiplicity_at(a::Real, energy::Real, lz::Real, q::Real, r::Real;
        atol::Real=_root_atol(_float_type(a, energy, lz, q, r)),
        rtol::Real=_root_rtol(_float_type(a, energy, lz, q, r)),
        energy_atol::Real=DEFAULT_ENERGY_ATOL,
        energy_rtol::Real=DEFAULT_ENERGY_RTOL)
    kwargs = (; energy_atol=energy_atol, energy_rtol=energy_rtol)
    coefficients = kerr_radial_coefficients(a, energy, lz, q; kwargs...)
    values = kerr_radial_derivatives(a, energy, lz, q, r; kwargs...)
    scales = _derivative_scales(coefficients, r)
    residuals = (values.R, values.R1, values.R2, values.R3, values.R4)
    near_zero(index) = abs(residuals[index]) <= atol + rtol * max(1.0, scales[index])

    near_zero(1) || return 0
    near_zero(2) || return 1
    near_zero(3) || return 2
    near_zero(4) || return 3
    near_zero(5) || return 4
    return 5
end

"""
    kerr_polar_theta_potential(a, E, Lz, Q, θ; axis_atol=1e-12)

``\\Theta_\\theta(\\theta) = Q - \\cos^2\\theta [a^2(1 - E^2) + L_z^2/\\sin^2\\theta]
= (d\\theta/d\\lambda)^2``. On the axis (``|\\sin\\theta| \\leq`` `axis_atol`) it
is finite only for ``L_z = 0``, where it equals ``Q - a^2(1 - E^2)``; otherwise it is `-Inf` there.
"""
function kerr_polar_theta_potential(a::Real, energy::Real, lz::Real, q::Real,
        theta::Real; axis_atol::Real=_tol(_float_type(a, energy, lz, q, theta), 1.0e-12))
    sine = sin(theta)
    cosine2 = cos(theta)^2
    if abs(sine) <= axis_atol
        _zero_lz(a, energy, lz, q) || return -_float_type(a, energy, lz, q, theta)(Inf)
        return q + cosine2 * a^2 * _e2m1(energy)
    end
    return q - cosine2 * (lz^2 / sine^2 - a^2 * _e2m1(energy))
end

"""
    kerr_polar_z_potential(a, E, Lz, Q, z)

``(dz/d\\lambda)^2`` for ``z = \\cos\\theta``:
``Q(1 - z^2) - z^2[L_z^2 + a^2(1 - E^2)(1 - z^2)]``, a polynomial in ``z`` that
stays regular on the axis.
"""
function kerr_polar_z_potential(a::Real, energy::Real, lz::Real, q::Real, z::Real)
    one_minus = 1 - z^2
    return q * one_minus - z^2 *
        (lz^2 - a^2 * _e2m1(energy) * one_minus)
end

"""
    kerr_polar_admissibility(a, E, Lz, Q; rtol=4eps(T))

Maximize ``(dz/d\\lambda)^2 = Q(1 - u) - L_z^2u + \\beta u(1 - u)``,
with ``u = \\cos^2\\theta`` and ``\\beta = a^2(E^2 - 1)``, over ``u \\in [0, 1]``
in closed form (regular at ``E = 1`` and ``a = 0``). The constants admit polar motion
(`admissible`) when the polar potential at one of the examined points is non-negative
to within what the constants themselves resolve:
``\\sum_j |\\partial\\Theta/\\partial c_j|\\,\\mathrm{ulp}(c_j)`` over
``(c_1,c_2,c_3,c_4) = (a,E,L_z,Q)``, plus `rtol` times
the size of the terms for the rounding of ``\\Theta`` (`T` the floating-point type of the
constants). A negative ``Q``, or a turning point short of the
axis, is therefore not admitted because its scale is small next to other terms. The axis
``u = 1`` (``\\Theta = -L_z^2``) is examined when ``L_z = 0`` and the motion reaches it. The result also carries
`max_value`, `max_cosine_squared`, the `tolerance` of that maximum and the examined
`candidates`.
"""
function kerr_polar_admissibility(a::Real, energy::Real, lz::Real, q::Real;
        rtol::Real=4eps(_float_type(a, energy, lz, q)))
    T = _float_type(a, energy, lz, q)
    beta = a^2 * _e2m1(energy)
    slope = beta - q - lz^2
    ulp(x) = iszero(x) ? zero(T) : eps(abs(T(x)))
    function candidate(u)
        w = u * (1 - u)
        value = q * (1 - u) - lz^2 * u + beta * w
        reach = abs(1 - u) * ulp(q) + 2abs(lz) * u * ulp(lz) +
            2abs(a) * abs(_e2m1(energy)) * w * ulp(a) + 2a^2 * abs(energy) * w * ulp(energy)
        terms = abs(q) * (1 - u) + lz^2 * u + abs(beta) * w
        return (u=T(u), value=T(value), tolerance=T(reach + rtol * terms))
    end
    candidates = [candidate(zero(T))]
    # the axis u = 1 (value −Lz²) counts when Lz = 0 and motion reaches it: Q + β ≥ 0
    if !_zero_lz(a, energy, lz, q) || q + beta >= -(ulp(q) + 2abs(a * _e2m1(energy)) * ulp(a) +
            2a^2 * abs(energy) * ulp(energy) + rtol * (abs(q) + abs(beta)))
        push!(candidates, candidate(one(T)))
    end
    if !iszero(beta)
        stationary = slope / (2 * beta)
        0 < stationary < 1 && push!(candidates, candidate(stationary))
    end
    best = candidates[argmax(getfield.(candidates, :value))]
    geometry = q < 0 && beta > 0 ? _polar_vortical_geometry(a, energy, lz, q) : nothing
    admissible = if geometry === nothing
        any(c -> c.value >= -c.tolerance, candidates)
    else
        # The reachable axis is an exact zero for Lz=0 even when the repeated
        # root's rounded center lies one ulp above u=1.
        (iszero(lz) && any(c -> c.u == 1.0 && iszero(c.value), candidates)) ||
            (0 < geometry.root_sum / 2 <= 1 &&
                (geometry.repeated || geometry.discriminant > 0))
    end
    return (
        admissible=admissible,
        max_value=best.value,
        max_cosine_squared=best.u,
        tolerance=best.tolerance,
        candidates=Tuple(candidates),
        reading=geometry === nothing ? :exact : geometry.reading,
    )
end

"""
    kerr_polar_sector_candidates(a, E, Lz, Q)

The polar sectors these constants allow, as a tuple of Symbols: `:pendular`, `:equatorial`,
`:equator_attractive`, `:vortical`, `:constant_latitude`, `:axis_crossing` or
`:axis_constant`. The sign of ``Q`` decides, and only ``Q = 0`` itself is equatorial. Where the
constants allow more than one (``Q = 0`` with ``L_z^2 < a^2(E^2 - 1)`` allows both `:equatorial` and
`:equator_attractive`), the keyword `polar_sector` of the constructors chooses the motion.
"""
function kerr_polar_sector_candidates(a::Real, energy::Real, lz::Real, q::Real)
    sectors = Symbol[]
    _on_axis(a, energy, lz, q) && return (:axis_constant,)
    # Lz = 0 with Q ≠ 0 and Q > a²(1-E²): the polar turning point is the pole itself, so
    # the motion passes over the axis. The pendular closed form degenerates there; the
    # vortical one (Q < 0) still describes the same motion and stays available on request.
    crossing = _axis_crossing(a, energy, lz, q)
    if crossing && !iszero(q)
        return q > 0 ? (:axis_crossing,) : (:axis_crossing, :vortical)
    end

    # the sign of Q decides; Q = 0 itself is the equatorial limit, where motion off the
    # equator (equator attractive) needs Lz² < a²(E² − 1)
    beta = a^2 * _e2m1(energy)
    if iszero(q)
        push!(sectors, :equatorial)
        lz^2 < beta && push!(sectors, :equator_attractive)
    elseif q > 0
        push!(sectors, :pendular)
    elseif beta > 0
        geometry = _polar_vortical_geometry(a, energy, lz, q)
        geometry.repeated || geometry.discriminant >= 0 || return (:unclassified_polar,)
        push!(sectors, geometry.repeated ? :constant_latitude : :vortical)
    end

    crossing && push!(sectors, :axis_crossing)
    isempty(sectors) && push!(sectors, :unclassified_polar)
    return Tuple(unique(sectors))
end

"""
The genuinely complex radial roots: the `complex_root_count` raw roots farthest from the
real axis. (A repeated real root comes out of the root finder as a nearly real pair
whose tiny imaginary parts must not be mistaken for a conjugate pair.)
"""
function _nonreal_roots(structure)
    n = structure.complex_root_count
    raw = sort(complex.(float.(collect(structure.raw_roots))); by=z -> -abs(imag(z)))
    return raw[1:min(n, length(raw))]
end

"""
    _double_root_factorization(a, E, Lz, Q, rc)

For E < 1 constants with an exterior double radial root near `rc`, polish `rc` to machine
precision (Newton on R' = 0, where the double root is simple) and return `(x1, rc, ra)`
with R(r) = (1 - E²)(r - x1)(r - rc)²(ra - r), the remaining roots taken from the
coefficients so that the factorization is exact.
"""
function _double_root_factorization(a, energy, lz, q, rc)
    coefficients = kerr_radial_coefficients(a, energy, lz, q)
    rc = _polish_root(coefficients, rc; order=1)
    c0, c1, c2, c3, c4 = coefficients
    kappa = -c4                                    # 1 - E²
    s = c3 / kappa - 2rc                           # x1 + ra
    p = -c0 / (kappa * rc^2)                       # x1 * ra
    disc = sqrt(max(s^2 - 4p, 0.0))
    return (x1=(s - disc) / 2, rc=rc, ra=(s + disc) / 2)
end

"""
    _polar_quadratic_roots(a, energy, lz, q)

Roots u = z^2 of Θ(z) = c u^2 - (q + lz^2 + c) u + q, with c = a^2(1 - E^2), in the
cancellation-free form (the textbook formula loses every digit of the small root when
|c| -> 0, i.e. a -> 0 or E -> 1). Returns the two roots and `c*u` for each (finite when
c = 0), so moduli and frequencies can be formed without dividing by c.
"""
function _polar_quadratic_roots(a, energy, lz, q)
    aa, ee, ll, qq = _wide.(float.((a, energy, lz, q)))
    beta = _wide_mul(_wide_mul(aa, aa), _wide_mul(_wide_sub(ee, _wide(1.0)), _wide_add(ee, _wide(1.0))))
    c = _wide_neg(beta)
    # s = Q + Lz² − β cancels when Lz² ≈ a²(E² − 1). The quadratic is solved in double-double
    # from the exact inputs and only the results are rounded, so identities of the exact
    # coefficients (Θ(1) = −Lz²) survive in the rounded roots.
    s = _wide_add(_wide_sub(_wide_mul(ll, ll), beta), qq)
    c1, s1, q1 = c[1] + c[2], s[1] + s[2], float(q)
    T = typeof(c1)
    if c1 < 0 < q1 || q1 < 0 < c1
        # disc = s² + 4|c q| has no cancellation; hypot keeps √disc when c q is subnormal
        sq = _wide(hypot(s1, 2 * sqrt(abs(c1)) * sqrt(abs(q1))))
        disc = s1^2 - 4 * c1 * q1
    else
        d = _wide_sub(_wide_mul(s, s), _wide_mul(_wide(4.0), _wide_mul(c, qq)))
        disc = d[1] + d[2]
        sq = disc > 0 ? _wide_sqrt(d) : _wide(zero(T))
    end
    big = _wide_mul(_wide(0.5), _wide_add(s, signbit(s1) ? _wide_neg(sq) : sq))   # |big| ≥ |s|/2
    big1 = big[1] + big[2]
    u_small = iszero(big1) ? zero(T) : q1 / big1
    u_big = iszero(c1) ? copysign(T(Inf), big1) : (u = _wide_div(big, c); u[1] + u[2])
    return (c=c1, disc=disc, u_small=u_small, u_big=u_big, cu_small=c1 * u_small, cu_big=big1)
end

# Coefficient form also serves interfaces that supply c directly.
_polar_quadratic_roots(c, lz, q) = _polar_quadratic_roots_from_sum(c, q + lz^2 + c, q)

function _polar_quadratic_roots_from_sum(c, s, q)
    T = _float_type(c, s, q)
    disc = s^2 - 4 * c * q
    # c and q of opposite signs: disc = s² + 4|c q| has no cancellation, and hypot keeps its
    # square root when s² or c q falls below the normal range (q down to the subnormals)
    sq = c < 0 < q || q < 0 < c ? hypot(s, 2 * sqrt(abs(c)) * sqrt(abs(q))) :
        sqrt(max(disc, 0.0))
    big = (s + copysign(sq, s)) / 2               # |big| >= |s|/2, no cancellation
    u_small = iszero(big) ? zero(T) : q / big
    u_big = iszero(c) ? copysign(T(Inf), big) : big / c
    return (c=c, disc=disc, u_small=u_small, u_big=u_big,
        cu_small=c * u_small, cu_big=big)
end

# The A reading uses the one-ulp input reach of D=(beta-Q-Lz^2)^2+4beta*Q,
# not a relative tolerance on its cancelling terms.
function _polar_vortical_geometry(a, energy, lz, q)
    roots = _polar_quadratic_roots(a, energy, lz, q)
    beta = -roots.c
    b = beta - q - lz^2
    ulp(x) = iszero(x) ? zero(b) : eps(abs(oftype(b, x)))
    dbeta = 2abs(a * _e2m1(energy)) * ulp(a) + 2a^2 * abs(energy) * ulp(energy)
    reach = abs(2b + 4q) * dbeta + abs(4beta - 2b) * ulp(q) +
        abs(4lz * b) * ulp(lz)
    exact = try
        isempty(_InputArithmetic.polar_discriminant(a, energy, lz, q))
    catch err
        err isa _InputArithmetic.Uncertified || rethrow()
        false
    end
    repeated = exact || abs(roots.disc) <= reach
    return (root_sum=b / beta, discriminant=roots.disc / beta^2,
        repeated=repeated, tolerance=reach / beta^2,
        reading=exact ? :exact : repeated ? :within_input_ulp : :exact)
end

function _polar_vortical_geometry(beta, lz, q; rtol=_constant_latitude_rtol(_float_type(beta, lz, q)))
    bb, ll, qq = _wide.(float.((beta, lz, q)))
    b = _wide_sub(_wide_sub(bb, qq), _wide_mul(ll, ll))
    d = _wide_add(_wide_mul(b, b), _wide_mul(_wide(4.0), _wide_mul(bb, qq)))
    slope, disc = b[1] + b[2], d[1] + d[2]
    ulp(x) = iszero(x) ? zero(slope) : eps(abs(oftype(slope, x)))
    reach = abs(2slope + 4q) * ulp(beta) + abs(4beta - 2slope) * ulp(q) +
        abs(4lz * slope) * ulp(lz)
    return (root_sum=slope / beta, discriminant=disc / beta^2,
        repeated=abs(disc) <= reach, tolerance=reach / beta^2,
        reading=iszero(disc) ? :uncertified_zero : abs(disc) <= reach ? :within_input_ulp : :exact)
end

# ---- tortoise coordinate ------------------------------------------------------------

"""
    kerr_rstar(a, r)

Tortoise coordinate with ``dr_*/dr = (r^2 + a^2)/\\Delta`` for ``r > r_+``:

```math
r_* = r + \\frac{2r_+}{d}\\ln\\frac{r-r_+}{2}
          - \\frac{2r_-}{d}\\ln\\frac{r-r_-}{2}, \\qquad d = r_+ - r_-.
```

At ``|a| = 1`` its limit is ``r + 2\\ln(r-1) - 2/(r-1) - 2\\ln 2``.
Returns `NaN` for ``r \\leq r_+``.
"""
function kerr_rstar(a::Real, r::Real)
    rp = _rplus(a)
    rm = _rminus(a)
    r <= rp && return NaN
    d = rp - rm
    u = r - 1
    if u > 200d
        # (2/d)[f(rp) - f(rm)] with f(x) = x log(r - x), expanded about x = 1: finite
        # as the horizons merge, r_* -> r + 2 log(r - 1) - 2/(r - 1) - 2 log 2 at |a| = 1.
        f1 = log(u) - 1 / u
        f3 = -1 / u^2 - 2r / u^3
        f5 = -6 / u^4 - 24r / u^5
        series = f1 + d^2 / 24 * f3 + d^4 / 1920 * f5
        series += _rstar_series_tail(r, u, d)
        return r + 2 * series - 2 * log(oftype(series, 2))
    end
    return r + 2 * rp / d * log((r - rp) / 2) -
           2 * rm / d * log((r - rm) / 2)
end

# The far-field expansions of r* and φ_H in d = r₊ − r₋ for u = r − 1 > 200d: odd powers k
# of d/2, each term at most (d/2u)² ≈ 2⁻¹⁷·³ of the previous one. The first three terms
# (k = 1, 3, 5) are written out above (Float64); `_far_terms(T)` odd terms reach eps(T), and
# the tails below add the terms k = 7, 9, … (none in Float64).
_far_terms(::Type{T}) where {T} = cld(precision(T) - 2, 17)
# r* (inside 2 × series): f⁽ᵏ⁾(1) (d/2)^(k−1)/k!, f(x) = x log(r − x), f⁽ᵏ⁾(1) = −(k−1)!/uᵏ − k(k−2)!/u^(k−1)
function _rstar_series_tail(r, u, d)
    T = _float_type(r, u, d)
    tail = zero(T)
    for j in 4:_far_terms(T)
        k = 2j - 1
        fk = -factorial(big(k - 1)) / T(u)^k - k * factorial(big(k - 2)) / T(u)^(k - 1)
        tail += T(fk * (T(d) / 2)^(k - 1) / factorial(big(k)))
    end
    return tail
end
# φ_H / a: −(d/2)^(k−1)/(k uᵏ)
function _azimuth_series_tail(u, d)
    T = _float_type(u, d)
    tail = zero(T)
    for j in 4:_far_terms(T)
        k = 2j - 1
        tail -= (T(d) / 2)^(k - 1) / (k * T(u)^k)
    end
    return tail
end

# Two-sided tortoise coordinate (inside r₊ as well), the horizon azimuth φ_H(r) with
# dφ_H/dr = a/Δ, and the horizon residues of the radial t and φ rates (P(r±) terms).
function _rstar_all(a, r)
    rp = _rplus(a)
    r > rp && return kerr_rstar(a, r)
    rm = _rminus(a)
    d = rp - rm
    inner = iszero(rm) ? zero(float(r)) : 2 * rm / d * log(abs(r - rm) / 2)
    return r + 2 * rp / d * log(abs(r - rp) / 2) - inner
end

function _horizon_azimuth(a, r)
    horizons = kerr_horizons(a)
    separation = horizons.rplus - horizons.rminus
    kerr_metric_limit(a) === :schwarzschild && return zero(_float_type(a, r))
    u = r - 1
    if abs(u) > 200separation
        # (a/d)[g(rp) - g(rm)] with g(x) = log|r - x|, expanded about x = 1; the extremal
        # limit is -a/(r - 1).
        series = -1 / u - separation^2 / 12 / u^3 - separation^4 / 80 / u^5
        series += _azimuth_series_tail(u, separation)
        return a * series
    end
    return a / separation * log(abs((r - horizons.rplus) /
        (r - horizons.rminus)))
end

function _radial_residues(a, energy, lz)
    horizons = kerr_horizons(a)
    separation = horizons.rplus - horizons.rminus
    abs(a) < 1 || error(
        "Horizon residues require r₊ > r₋ (|a| < 1); at |a| = 1 the two horizons coincide.")
    pplus = kerr_radial_momentum(a, energy, lz, horizons.rplus)
    pminus = kerr_radial_momentum(a, energy, lz, horizons.rminus)
    return (
        rplus=horizons.rplus,
        rminus=horizons.rminus,
        pplus=pplus,
        pminus=pminus,
        c_phi_plus=a * pplus / separation,
        c_phi_minus=-a * pminus / separation,
        c_t_plus=2 * horizons.rplus * pplus / separation,
        c_t_minus=-2 * horizons.rminus * pminus / separation,
    )
end
