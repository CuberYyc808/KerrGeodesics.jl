# Kerr metric functions shared by everything: horizons, energy regime, radial and polar
# potentials and their roots, and the tortoise coordinate.

const DEFAULT_CLASSIFICATION_ATOL = 64 * eps(Float64)
const DEFAULT_ENERGY_ATOL = 0.0
const DEFAULT_ENERGY_RTOL = 0.0

"""
    kerr_delta(a, r)

Δ(r) = r² − 2r + a² = (r − r₊)(r − r₋), evaluated in the product form, which keeps its
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

Horizon radii `(rplus, rminus)` = 1 ± √(1 − a²) for |a| ≤ 1 (G = c = M = 1).
"""
function kerr_horizons(a::Real; atol::Real=DEFAULT_CLASSIFICATION_ATOL)
    abs(a) <= 1 + atol || throw(DomainError(a, "Kerr spin must satisfy |a|<=1."))
    spin = clamp(float(a), -1.0, 1.0)
    return (rplus=_rplus(spin), rminus=_rminus(spin))
end

"""
    kerr_metric_limit(a; atol=64eps(), near_extremal_threshold=1e-6)

`:schwarzschild` (|a| ≤ `atol`), `:extremal` (||a| − 1| ≤ `atol`), `:near_extremal`
(1 − |a| ≤ `near_extremal_threshold`) or `:subextremal`, tested in that order.
"""
function kerr_metric_limit(a::Real;
        atol::Real=DEFAULT_CLASSIFICATION_ATOL,
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

`:elliptic`, `:parabolic` or `:hyperbolic` by the sign of E² − 1, the r⁴ coefficient of R
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

`+1` for E ≥ 0 and `-1` for E < 0. Future-directed motion with E < 0 exists only inside the
ergoregion: the Trapped class.
"""
kerr_energy_sign(energy::Real) = energy < 0 ? -1 : 1

"""
    kerr_radial_momentum(a, E, Lz, r)

P(r) = E(r² + a²) − a Lz. It enters the radial potential, R = P² − Δ[r² + (Lz − aE)² + Q],
and the rates dt/dλ and dφ/dλ; on the horizon R(r₊) = P(r₊)².
"""
kerr_radial_momentum(a::Real, energy::Real, lz::Real, r::Real) =
    energy * (r - 1) * (r + 1) + (energy * (1 + a^2) - a * lz)   # r² − 1 accurate near r = 1

"""
    kerr_axis_carter_q(a, E)

Carter constant Q = a²(1 − E²) of a timelike geodesic along the spin axis (Lz = 0).
"""
kerr_axis_carter_q(a::Real, energy::Real) = -a^2 * _e2m1(energy)

"""
    kerr_axis_radial_potential(a, E, r)

Radial potential on the spin axis (Lz = 0, Q = a²(1 − E²)) in factored form,
R = (r² + a²)[E²(r² + a²) − Δ].
"""
function kerr_axis_radial_potential(a::Real, energy::Real, r::Real)
    sigma = r^2 + a^2
    return sigma * (energy^2 * sigma - kerr_delta(a, r))
end

"""
    kerr_radial_coefficients(a, E, Lz, Q; energy_atol=0, energy_rtol=0)

Coefficients (c₀, c₁, c₂, c₃, c₄) of R(r) = Σ cₖ rᵏ, in ascending order. c₄ = E² − 1 is
set to zero when `kerr_energy_regime(E; atol=energy_atol, rtol=energy_rtol)` is
`:parabolic`, so R is built as a cubic; with the default zero tolerances this is E² = 1
exactly.
"""
function kerr_radial_coefficients(a::Real, energy::Real, lz::Real, q::Real;
        energy_atol::Real=DEFAULT_ENERGY_ATOL,
        energy_rtol::Real=DEFAULT_ENERGY_RTOL)
    regime = kerr_energy_regime(energy; atol=energy_atol, rtol=energy_rtol)
    c4 = regime === :parabolic ? 0.0 : _e2m1(energy)
    return (
        -a^2 * q,
        2 * (a * energy - lz)^2 + 2q,
        -(q + lz^2 - a^2 * _e2m1(energy)),
        2.0,
        c4,
    )
end

"""
    kerr_radial_polynomial(a, E, Lz, Q; energy_atol=0, energy_rtol=0)

R(r) as a `Polynomial` (Polynomials.jl) built from `kerr_radial_coefficients`; it is a
cubic when E² = 1.
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

Radial potential of a timelike Kerr geodesic with Carter constant `Q`,
R = [E(r² + a²) − aLz]² − Δ[r² + (Lz − aE)² + Q]. It is summed from its coefficients
(`kerr_radial_coefficients`), which keeps full relative precision at large r, where the
factored form would cancel two terms of size E²r⁴.
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
    kerr_root_multiplicity_at(a, E, Lz, Q, r; atol=1e-9, rtol=1e-9, energy_atol=0,
                              energy_rtol=0)

Multiplicity of the radius `r` as a root of R: the number of consecutive values R, R′, R″,
R‴, R⁗ at `r`, starting from R, that vanish within `atol + rtol·max(1, s)`, where `s` is the
sum of the magnitudes of that derivative's terms (0 when R(r) ≠ 0).
"""
function kerr_root_multiplicity_at(a::Real, energy::Real, lz::Real, q::Real, r::Real;
        atol::Real=1.0e-9,
        rtol::Real=1.0e-9,
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

Θ(θ) = Q − cos²θ [a²(1 − E²) + Lz²/sin²θ] = (dθ/dλ)². On the axis (sin θ = 0) it is finite
only for Lz = 0, where it equals Q − a²(1 − E²); otherwise it is `-Inf` there.
"""
function kerr_polar_theta_potential(a::Real, energy::Real, lz::Real, q::Real,
        theta::Real; axis_atol::Real=1.0e-12)
    sine = sin(theta)
    cosine2 = cos(theta)^2
    if abs(sine) <= axis_atol
        abs(lz) <= axis_atol || return -Inf
        return q + cosine2 * a^2 * _e2m1(energy)
    end
    return q - cosine2 * (lz^2 / sine^2 - a^2 * _e2m1(energy))
end

"""
    kerr_polar_z_potential(a, E, Lz, Q, z)

(dz/dλ)² for z = cos θ: Q(1 − z²) − z²[Lz² + a²(1 − E²)(1 − z²)], a polynomial in z that
stays regular on the axis.
"""
function kerr_polar_z_potential(a::Real, energy::Real, lz::Real, q::Real, z::Real)
    one_minus = 1 - z^2
    return q * one_minus - z^2 *
        (lz^2 - a^2 * _e2m1(energy) * one_minus)
end

"""
    kerr_polar_admissibility(a, E, Lz, Q; atol=1e-12, rtol=1e-10)

Maximize (dz/dλ)² = Q + (β − Q − Lz²)u − βu², u = cos²θ, β = a²(E² − 1), over u ∈ [0, 1] in
closed form (regular at E = 1 and a = 0). The constants admit polar motion (`admissible`)
when the maximum is at least −`tolerance`; the result also carries `max_value`,
`max_cosine_squared`, `tolerance` and the examined `candidates`.
"""
function kerr_polar_admissibility(a::Real, energy::Real, lz::Real, q::Real;
        atol::Real=1.0e-12,
        rtol::Real=1.0e-10)
    beta = a^2 * _e2m1(energy)
    value(u) = q + (beta - q - lz^2) * u - beta * u^2
    tolerance = atol + rtol * max(1.0, abs(q), abs(beta), lz^2)
    candidates = [(u=0.0, value=value(0.0))]
    axis_theta = q + beta
    if abs(lz) > atol || axis_theta >= -tolerance
        push!(candidates, (u=1.0, value=value(1.0)))
    end
    if abs(beta) > atol
        stationary = (beta - q - lz^2) / (2 * beta)
        if 0 < stationary < 1
            push!(candidates, (u=float(stationary), value=value(stationary)))
        end
    end
    best = candidates[argmax(getfield.(candidates, :value))]
    return (
        admissible=best.value >= -tolerance,
        max_value=best.value,
        max_cosine_squared=best.u,
        tolerance=tolerance,
        candidates=Tuple(candidates),
    )
end

# |Lz| below which the orbit is treated as passing over the axis: with Lz^2 < eps*Q the
# polar turning point z_+ = 1 - O(Lz^2/Q) rounds to the pole, so the generic pendular
# forms have no digits left and the Lz = 0 (axis-crossing) forms are the exact ones.
# (z_+ = 1 - O(Lz^2/max(Q, |c|)), c = a^2(1 - E^2)).
_axis_lz_tolerance(q, c=0.0; atol=1.0e-12) =
    max(atol, 2 * sqrt(eps(1.0)) * sqrt(max(abs(q), abs(c))))

"""
    kerr_polar_sector_candidates(a, E, Lz, Q; atol=1e-12)

The polar sectors these constants allow, as a tuple of Symbols: `:pendular`, `:equatorial`,
`:equator_attractive`, `:vortical`, `:constant_latitude`, `:axis_crossing` or
`:axis_constant`. Where the constants allow more than one (Q = 0 with E > 1 allows both
`:equatorial` and `:equator_attractive`), the keyword `polar_sector` of the constructors
chooses the motion.
"""
function kerr_polar_sector_candidates(a::Real, energy::Real, lz::Real, q::Real;
        atol::Real=1.0e-12)
    sectors = Symbol[]
    lz_axis = _axis_lz_tolerance(q, -a^2 * _e2m1(energy); atol=atol)
    axis_value = q + a^2 * _e2m1(energy)
    axis_tolerance = atol * max(1.0, abs(q), abs(a^2 * _e2m1(energy)))
    if abs(lz) <= lz_axis && abs(axis_value) <= axis_tolerance
        return (:axis_constant,)
    end
    # Lz = 0 with Q ≠ 0 and Q > a²(1-E²): the polar turning point is the pole itself, so
    # the motion passes over the axis. The pendular closed form degenerates there; the
    # vortical one (Q < 0) still describes the same motion and stays available on request.
    if abs(lz) <= lz_axis && axis_value > axis_tolerance && abs(q) > atol
        return q > 0 ? (:axis_crossing,) : (:axis_crossing, :vortical)
    end

    if abs(q) <= atol
        push!(sectors, :equatorial)
        energy^2 > 1 + atol && push!(sectors, :equator_attractive)
    elseif q > 0
        push!(sectors, :pendular)
    elseif energy^2 > 1 + atol
        beta = a^2 * _e2m1(energy)
        if beta > atol
            root_sum = (beta - q - lz^2) / beta
            root_product = -q / beta
            discriminant = root_sum^2 - 4 * root_product
            tolerance = 128 * atol * max(
                1.0, root_sum^2, abs(root_product))
            push!(sectors, abs(discriminant) <= tolerance ?
                :constant_latitude : :vortical)
        else
            push!(sectors, :unclassified_polar)
        end
    end

    if abs(lz) <= lz_axis && axis_value > axis_tolerance
        push!(sectors, :axis_crossing)
    end
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
    raw = sort(ComplexF64.(collect(structure.raw_roots)); by=z -> -abs(imag(z)))
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
    for _ in 1:8
        d = kerr_radial_derivatives(a, energy, lz, q, rc)
        step = d.R1 / d.R2
        rc -= step
        abs(step) <= 4eps(rc) * max(1.0, abs(rc)) && break
    end
    c0, c1, c2, c3, c4 = kerr_radial_coefficients(a, energy, lz, q)
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
    c = -a^2 * _e2m1(energy)
    s = q + lz^2 + c
    disc = s^2 - 4 * c * q
    sq = sqrt(max(disc, 0.0))
    big = (s + copysign(sq, s)) / 2               # |big| >= |s|/2, no cancellation
    u_small = iszero(big) ? 0.0 : q / big
    u_big = iszero(c) ? copysign(Inf, big) : big / c
    return (c=c, disc=disc, u_small=u_small, u_big=u_big,
        cu_small=c * u_small, cu_big=big)
end

"""
    _elliptic_D(φ, m)

D(φ|m) = ∫_0^φ sin²θ / sqrt(1 - m sin²θ) dθ = (F(φ|m) - E(φ|m))/m, evaluated without the
cancellation of F - E at small m (a series in m for |m| < 0.05). Polar t pieces are
proportional to (F - E)/(1 - E^2) and would otherwise lose ~eps/(1 - E^2) as E -> 1.
"""
function _elliptic_D(φ, m)
    if abs(m) >= 0.05
        return (Elliptic.F(φ, m) - Elliptic.E(φ, m)) / m
    end
    s, c = sincos(φ)
    S = φ                                     # S_j = ∫ sin^{2j}
    b = 1.0                                   # binomial(2n,n)/4^n
    total = 0.0
    mn = 1.0
    for n in 0:40
        j = n + 1
        S = ((2j - 1) * S - s^(2j - 1) * c) / (2j)
        term = b * mn * S
        total += term
        abs(term) <= 1.0e-17 * abs(total) && n > 2 && break
        b *= (2n + 1) / (2n + 2)
        mn *= m
    end
    return total
end

# ---- tortoise coordinate ------------------------------------------------------------

"""
    kerr_rstar(a, r)

Tortoise coordinate with dr*/dr = (r² + a²)/Δ for r > r₊,
r* = r + (2r₊/d) log((r − r₊)/2) − (2r₋/d) log((r − r₋)/2),  d = r₊ − r₋,
with the limit r + 2 log(r − 1) − 2/(r − 1) − 2 log 2 at |a| = 1. Returns `NaN` for r ≤ r₊.
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
        return r + 2 * (f1 + d^2 / 24 * f3 + d^4 / 1920 * f5) - 2 * log(2)
    end
    return r + 2 * rp / d * log((r - rp) / 2) -
           2 * rm / d * log((r - rm) / 2)
end

# Two-sided tortoise coordinate (inside r₊ as well), the horizon azimuth φ_H(r) with
# dφ_H/dr = a/Δ, and the horizon residues of the radial t and φ rates (P(r±) terms).
function _rstar_all(a, r)
    rp = _rplus(a)
    r > rp && return kerr_rstar(a, r)
    rm = _rminus(a)
    d = rp - rm
    inner = iszero(rm) ? 0.0 : 2 * rm / d * log(abs(r - rm) / 2)
    return r + 2 * rp / d * log(abs(r - rp) / 2) - inner
end

function _horizon_azimuth(a, r)
    horizons = kerr_horizons(a)
    separation = horizons.rplus - horizons.rminus
    abs(a) <= 1.0e-15 && return 0.0
    u = r - 1
    if abs(u) > 200separation
        # (a/d)[g(rp) - g(rm)] with g(x) = log|r - x|, expanded about x = 1; the extremal
        # limit is -a/(r - 1).
        return a * (-1 / u - separation^2 / 12 / u^3 - separation^4 / 80 / u^5)
    end
    return a / separation * log(abs((r - horizons.rplus) /
        (r - horizons.rminus)))
end

function _radial_residues(a, energy, lz)
    horizons = kerr_horizons(a)
    separation = horizons.rplus - horizons.rminus
    separation > 1.0e-12 || error(
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
        c_t_plus=2.0 * horizons.rplus * pplus / separation,
        c_t_minus=-2.0 * horizons.rminus * pminus / separation,
    )
end
