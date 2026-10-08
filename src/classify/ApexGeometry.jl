# Radial root geometry of APEX input. (p, e, x) fixes the turning points r1 = p/(1 − e) and
# r2 = p/(1 + e) exactly; the constants (E, Lz, Q) computed from them are rounded. Near the
# separatrix r2 − r3 is small and a root structure recomputed from the rounded constants can
# merge r2 and r3 into a repeated root, so the bound orbit would be lost. The geometry below
# keeps r1 and r2 and divides R(r) of the same constants, in double-double, down to its inner quadratic factor; the
# Stable member reads its roots from it. Every other member, and constants input, keep the
# root structure of the constants.

# Divide c (ascending coefficients) by (r − ρ) from the constant term up: stable when ρ is
# the largest root in magnitude.
function _deflate_largest(c::AbstractVector{<:Real}, ρ)
    n = length(c) - 1
    q = zeros(float(eltype(c)), n)
    q[1] = -c[1] / ρ
    for k in 2:n
        q[k] = (q[k - 1] - c[k]) / ρ
    end
    return q
end

# The inner quadratic c2 r² + c1 r + c0 of R(r)/((r − r1)(r − r2)) for the rounded constants,
# with the coefficients of R formed and divided in double-double: near the separatrix the
# quadratic's roots are needed to a small fraction of r2 − r3, which a working-precision
# division leaves to the cancellation in its steps. Returns its coefficients (rounded) and the
# roots r3 ≥ r4.
function _apex_inner_quadratic(a, energy, lz, q, r1, r2)
    c = _wide_radial_coefficients(a, energy, lz, q)
    divide(c, ρ) = begin
        out = Vector{typeof(c[1])}(undef, length(c) - 1)
        out[1] = _wide_div(_wide_neg(c[1]), _wide(ρ))
        for k in 2:length(out)
            out[k] = _wide_div(_wide_sub(out[k - 1], c[k]), _wide(ρ))
        end
        out
    end
    c0, c1, c2 = divide(divide(collect(c), r1), r2)
    disc = _wide_sub(_wide_mul(c1, c1), _wide_mul(_wide_mul(_wide(4one(r1)), c2), c0))
    rounded(v) = v[1] + v[2]
    coefficients = (rounded(c0), rounded(c1), rounded(c2), rounded(disc))
    disc[1] > 0 || return coefficients, (oftype(r1, NaN), oftype(r1, NaN))
    root = _wide_sqrt(disc)
    big = _wide_mul(_wide(-one(r1) / 2), c1[1] >= 0 ? _wide_add(c1, root) : _wide_sub(c1, root))
    ra, rb = rounded(_wide_div(big, c2)), rounded(_wide_div(c0, big))
    return coefficients, (max(ra, rb), min(ra, rb))
end

"""
    _apex_root_geometry(a, p, e, x, E, Lz, Q)

The radial roots of a bound eccentric APEX orbit with four separated simple roots, one inside
the outer horizon: `(accepted, reason, structure, diagnostics)`. `structure` (a root
structure as `kerr_geo_root_structure` returns it, with the roots r1 = p/(1 − e),
r2 = p/(1 + e) and r3 > r4 of the deflated quadratic) is given when `accepted`; otherwise
`reason` names the first condition that fails and the orbit keeps the root structure of the
constants. The conditions are working-precision checks:

- `:outside_domain`: |a| < 1, 0 < e < 1, 0 < E < 1, Q ≥ 0 do not hold;
- `:horizon_root`: P(r₊) = 0 (the horizon is a root of R);
- `:turning_points_not_exterior`: r1 > r2 > r₊ does not hold;
- `:discriminant_unresolved`: the quadratic factor's discriminant is not positive beyond the
  rounding of its two terms;
- `:roots_not_ordered`: r2 > r3 > r4 does not hold;
- `:gap_unresolved`: r2 − r3 is within the rounding of r3 and r2;
- `:backward_contract`: R of the constants and the factored model differ by more than
  256 eps (root residuals or coefficients, relative);
- `:topology`: not one root inside and three outside the outer horizon.
"""
function _apex_root_geometry(a, p, e, x, energy, lz, q)
    T = _float_type(a, p, e, x)
    reject(reason, diagnostics=NamedTuple()) =
        (accepted=false, reason=reason, structure=nothing, diagnostics=diagnostics)
    abs(a) < 1 && 0 < e < 1 && 0 < energy < 1 && q >= 0 || return reject(:outside_domain)
    _horizon_root(a, energy, lz) && return reject(:horizon_root)
    r1, r2 = p / (1 - e), p / (1 + e)
    rplus = _rplus(a)
    r1 > r2 > rplus || return reject(:turning_points_not_exterior)
    c = collect(T, kerr_radial_coefficients(a, energy, lz, q))
    (c0, c1, c2, disc), (r3, r4) = _apex_inner_quadratic(a, energy, lz, q, r1, r2)
    unit = eps(T)
    disc_bound = 32unit * (abs(c1^2) + abs(4c2 * c0))
    disc > disc_bound || return reject(:discriminant_unresolved, (disc=disc, disc_bound=disc_bound))
    all(isfinite, (r3, r4)) && r2 > r3 > r4 || return reject(:roots_not_ordered)
    # sensitivity of r3 to the rounding of the quadratic's coefficients, plus the spacing of
    # the floating-point numbers at r3 and r2
    gap_bound = 64unit * (abs(c0) + abs(c1 * r3) + abs(c2 * r3^2)) / abs(c1 + 2c2 * r3) +
        32eps(r3) + 32eps(r2)
    r2 - r3 > gap_bound || return reject(:gap_unresolved, (gap=r2 - r3, gap_bound=gap_bound))
    roots = (r1, r2, r3, r4)
    evaluator = _radial_root_evaluator(a, energy, lz, q)
    scaled(r) = abs(first(evaluator(r))) / max(one(T), sum(abs(c[j + 1]) * abs(r)^j for j in 0:4))
    backward = maximum(scaled, roots)
    # coefficients of c4 ∏(r − rᵢ) against those of the constants
    model = T[c[end]]
    for r in roots
        next = zeros(T, length(model) + 1)
        for j in eachindex(model)
            next[j] -= r * model[j]
            next[j + 1] += model[j]
        end
        model = next
    end
    mismatch = maximum(abs(model[j] - c[j]) / max(one(T), abs(c[j])) for j in eachindex(c))
    diagnostics = (gap=r2 - r3, gap_bound=gap_bound, disc=disc, disc_bound=disc_bound,
        root_backward=backward, coefficient_backward=mismatch, bound=256unit)
    max(backward, mismatch) <= 256unit || return reject(:backward_contract, diagnostics)
    items = map(reverse(roots)) do r
        v = kerr_radial_derivatives(a, energy, lz, q, r)
        (radius=r, multiplicity=1, source=:apex_turning_points, reading=:apex_turning_points,
            residuals=(R=v.R, R1=v.R1, R2=v.R2, R3=v.R3))
    end
    below = Tuple(item for item in items if item.radius < rplus)
    exterior = Tuple(item for item in items if item.radius > rplus)
    length(below) == 1 && length(exterior) == 3 || return reject(:topology, diagnostics)
    structure = (polynomial=kerr_radial_polynomial(a, energy, lz, q), degree=4,
        raw_roots=collect(Complex{T}, reverse(roots)), real_roots=items,
        below_horizon=below, horizon_coincident=(), exterior=exterior,
        complex_root_count=0, multiplicity_sum=4,
        root_tolerance=(atol=_root_atol(T), rtol=_root_rtol(T)), near_repeated=(),
        apex_turning_points=(roots=roots, diagnostics...))
    return (accepted=true, reason=:accepted, structure=structure, diagnostics=diagnostics)
end
