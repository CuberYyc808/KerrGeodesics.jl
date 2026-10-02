# Degeneracy criteria: the named tolerances and predicates that decide when constants sit on a
# boundary between cases (E = 1, Q = 0, Lz = 0 and the spin axis, a radial root on the outer
# horizon, the spin limits a = 0 and |a| = 1, a constant latitude) and when a Mino time or a
# radius sits on a domain endpoint. Every site decides through these, so each criterion has
# one definition. Q = 0 is exact (`iszero(q)`): the sign of Q selects the polar sector. E = 1
# is exact as well (`kerr_energy_regime`, core/Metric.jl): the sign of E² − 1 selects the
# radial formulas.

# Spin limits: a = 0 (Schwarzschild) and |a| = 1 (extremal) within a few ulps
# (`kerr_metric_limit`). APEX input alone rounds a spin within SPIN_SNAP_TOL of ±1 to ±1.
const DEFAULT_CLASSIFICATION_ATOL = 64 * eps(Float64)
const SPIN_SNAP_TOL = 8 * eps(Float64)

# E = 1 is exact by default; the keywords of `kerr_energy_regime` and of the radial
# coefficients remain for callers that snap the quartic coefficient of R on purpose.
const DEFAULT_ENERGY_ATOL = 0.0
const DEFAULT_ENERGY_RTOL = 0.0

# Radial roots: two roots this close are one repeated root, a root this close to r₊ lies on
# the horizon, and a probe radius this close to a root sits on it (`kerr_geo_root_structure`,
# `kerr_root_multiplicity_at`).
const ROOT_ATOL = 1.0e-12
const ROOT_RTOL = 1.0e-12

# A root of multiplicity m at r: the Float64 constants must lie within one ulp per component
# of the manifold where R, R′, …, R^(m−1) vanish at r. To first order that is
#     |R^(k)(r)| ≤ Σ_θ |∂R^(k)/∂θ (r)| ulp(θ),   θ = (a, E, Lz, Q),   k = 0, …, m − 1,
# except k = `located`, the derivative whose zero placed r (it vanishes up to the rounding of
# r, which is not a property of the constants);
# with ∂R^(k)/∂θ from the coefficients cⱼ(θ) of R; an exactly zero constant carries no
# rounding. Distinct roots of the given constants (a resolved gap, a complex pair, a simple
# root next to a double one) are therefore not merged, while constants rounded from a point
# of the manifold are. R and R′ are evaluated in double-double, so their own rounding does
# not enter.
function _repeated_radial_zero(a,E,L,Q,r; multiplicity=2, located=multiplicity-1)
    value, slope = _radial_root_evaluator(a,E,L,Q)(r)
    higher = kerr_radial_derivatives(a,E,L,Q,r)
    values = (value, slope, higher.R2, higher.R3, higher.R4)
    for k in 0:min(multiplicity-1, 4)
        k == located && continue
        abs(values[k+1]) <= _radial_input_reach(a,E,L,Q,r,k) || return false
    end
    return true
end

function _radial_input_reach(a,E,L,Q,r,k)
    w = a*E-L
    # ∂(c₀, …, c₄)/∂θ for c = (−a²Q, 2(aE − Lz)² + 2Q, −(Q + Lz² − a²(E² − 1)), 2, E² − 1)
    partials = ((-2a*Q, 4E*w, 2a*_e2m1(E), 0.0, 0.0), (0.0, 4a*w, 2a^2*E, 0.0, 2E),
                (0.0, -4w, -2L, 0.0, 0.0), (-a^2, 2.0, -1.0, 0.0, 0.0))
    ulps = map(x -> iszero(x) ? 0.0 : eps(abs(float(x))), (a, E, L, Q))
    reach = 0.0
    for (dc, u) in zip(partials, ulps)
        reach += abs(sum(dc[j+1]*factorial(j)/factorial(j-k)*r^(j-k) for j in k:4)) * u
    end
    return reach
end

# A conservative exact-zero certificate at the represented radius, after locating
# R^(m-1)=0. An uncertifiable radius may still have the explicit A reading.
function _exact_repeated_radial_zero(a,E,L,Q,r,m)
    c = _wide_radial_coefficients(a,E,L,Q)
    coefficients = ntuple(k -> _wide_derivative_coefficients(c,k), 4)
    value(k, x) = k == 0 ? _wide_evalpoly(x, c) : _wide_evalpoly(x, coefficients[k])
    for _ in 1:8
        step = value(m-1, r) / value(m, r)
        isfinite(step) || break
        r -= step
        abs(step) <= eps(r) && break
    end
    return _InputArithmetic.repeated_at(a,E,L,Q,r,m), r
end

# P(r₊) = 0: r₊ itself is a root of R (the H tier and, at |a| = 1, the X tier). P(r₊) is a sum
# of two products and is exact to a few ulps, so a few ulps is the test; a root numerically at
# r₊ with P(r₊) ≠ 0 bounds the allowed sliver above the horizon instead (`kerr_geo_root_structure`).
const HORIZON_MOMENTUM_TOL = DEFAULT_CLASSIFICATION_ATOL
_horizon_root(a, energy, lz) = abs(kerr_radial_momentum(a, energy, lz, _rplus(a))) <=
    HORIZON_MOMENTUM_TOL * max(1.0, abs(energy), abs(lz))

# Lz = 0: |Lz| below which the polar turning point rounds to the pole, z₊ = 1 − O(Lz²/max(Q, |c|))
# with c = a²(1 − E²), so the generic pendular forms have no digits left and the Lz = 0
# (axis-crossing) forms are the exact ones. The bound scales with Q and c, so a small Q does
# not turn a small Lz into an axis orbit.
_axis_lz_tolerance(q, c) = 2 * sqrt(eps(1.0)) * sqrt(max(abs(q), abs(c)))
_zero_lz(a, energy, lz, q) = abs(lz) <= _axis_lz_tolerance(q, kerr_axis_carter_q(a, energy))

# Motion along the spin axis: Lz = 0 and Q = a²(1 − E²).
function _on_axis(a, energy, lz, q)
    c = kerr_axis_carter_q(a, energy)
    return abs(lz) <= _axis_lz_tolerance(q, c) && abs(q - c) <= ROOT_RTOL * max(abs(q), abs(c))
end

# Motion over the axis: Lz = 0 with Q above a²(1 − E²), so the polar turning point is the pole.
function _axis_crossing(a, energy, lz, q)
    c = kerr_axis_carter_q(a, energy)
    return abs(lz) <= _axis_lz_tolerance(q, c) && q - c > ROOT_RTOL * max(abs(q), abs(c))
end

# Constant latitude: the two roots of the vortical polar polynomial coincide, relative to the
# two terms of its discriminant (`_polar_vortical_geometry`).
const CONSTANT_LATITUDE_RTOL = 128 * ROOT_RTOL

# A computed polar turning root may land this far above 1: Θ(1) = −Lz² puts z₊ ≤ 1, and with
# Lz² ≪ Q the rounded root can sit an ulp above 1.
const POLAR_ROOT_SLACK = 4 * eps(Float64)

# A Mino time this close outside a closed domain endpoint is that endpoint, and a radius this
# close to a turning point sits on it (input rounding).
const MINO_ENDPOINT_TOL = 2.0e-12
const RADIUS_TOL = 2.0e-11

# The horizon-root quadratic h(z) = h2 z² + h1 z + h0 of |a| = 1 with P(r₊) = 0 (h2 = 3E² − 1 − Q,
# h1 = 4E² − 2, h0 = E² − 1): a coefficient this small counts as zero when the family is
# classified (h linear, the triple horizon root; both h1 and h2 zero, the excluded quadruple
# root). The radial model itself is uniform in h2 and h0 and needs no such threshold.
const HORIZON_ROOT_COEFFICIENT_TOL = 1.0e-11
