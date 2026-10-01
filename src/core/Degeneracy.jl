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

# A small residual of the expanded polynomial need not be a repeated root:
# near the extremal horizon its O(1) terms can hide a nonzero P_H^2.
# Require the metric form to vanish within its arithmetic/input rounding scale.
function _repeated_radial_zero(a,E,L,Q,r)
    P=kerr_radial_momentum(a,E,L,r)
    D=(r-1)^2+(a-1)*(a+1)
    K=r^2+(L-a*E)^2+Q
    value=muladd(P,P,-D*K)
    pscale=abs(E)*(r^2+a^2)+abs(a*L)
    dscale=abs((r-1)^2)+abs((a-1)*(a+1))
    kscale=r^2+(abs(L)+abs(a*E))^2+abs(Q)
    scale=2abs(P)*pscale+P^2+dscale*kscale+abs(D*K)
    # Eight rounded operations bound the longest path through P, Delta and K.
    # A classification-wide tolerance would merge resolved nearly circular roots.
    gamma=8eps(Float64)/(1-8eps(Float64))
    return abs(value)<=gamma*scale
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
