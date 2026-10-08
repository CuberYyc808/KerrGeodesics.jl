# Plunge reference API: Mino time as a function of radius (`lambda_of_r`) for each root class.

# Every Mino time is formed from the distances of the radius to the two roots that bracket it;
# λ_H takes them for r₊, `gaps = (outer − r₊, r₊ − inner)`, from the exact inputs
# (`_horizon_gap`): a root next to r₊ (P(r₊) → 0) leaves no digits in a rounded difference.

# F(φ|m) for φ ∈ [0, π/2] from s = sin²φ, c = cos²φ and m1 = 1 − m, in Carlson's form
# √s R_F(c, c + m1 s, 1): no angle and no K are formed, so φ → π/2 and m → 1 lose nothing.
_legendre_f(s, c, m1) = sqrt(s) * _carlson_rf(c, c + m1 * s, one(m1))

struct _HorizonMinoMap{T,F} <: Function
    radial::F
    horizon::T
    lambda_horizon::T
end
_HorizonMinoMap(radial::F, horizon, lambda_horizon) where {F} =
    _HorizonMinoMap{_float_type(horizon, lambda_horizon),F}(radial, horizon, lambda_horizon)
(map::_HorizonMinoMap)(r) = r == map.horizon ? map.lambda_horizon : map.radial(r)

function lambda_of_r_real(a, E, roots; gaps=(roots[2] - _rplus(a), _rplus(a) - roots[1]))
    r4, r3, r2, r1 = roots
    m1 = (r2 - r3) * (r1 - r4) / ((r1 - r3) * (r2 - r4))
    norm = sqrt(-_e2m1(E) * (r1 - r3) * (r2 - r4))
    # from r3 down to r3 − d3 = r4 + d4: sin²φ = d3 (r2 − r4)/((r2 − r)(r3 − r4)), cos²φ likewise
    function λ(d3, d4)
        den = (r2 - (r3 - d3)) * (r3 - r4)
        return 2 * _legendre_f(d3 * (r2 - r4) / den, d4 * (r2 - r3) / den, m1) / norm
    end
    function λ_of_r(r)
        r4 <= r <= r3 && return λ(r3 - r, r - r4)
        @info("r = $r is out of the plunge region between r4 = $r4 and r3 = $r3.")
    end
    Λr_max = λ_of_r(r4)
    d3, d4 = gaps
    (d3 >= 0 && d4 >= 0) || return Λr_max, λ_of_r(_rplus(a)), λ_of_r   # r₊ outside the range
    λH = λ(d3, d4)
    return Λr_max, λH, _HorizonMinoMap(λ_of_r, _rplus(a), λH)
end

function lambda_of_r_complex(a, E, roots; gaps=(roots[1] - _rplus(a), _rplus(a) - roots[2]),
        pair)
    r1, r2, A, B = roots
    ρ, η = pair
    # X − x for X = √(x² + η²), without cancellation; m1 = ((A + B)² − (r1 − r2)²)/(4AB)
    excess(X, x) = x > 0 ? η^2 / (X + x) : X - x
    m1 = (excess(A, r1 - ρ) + excess(B, ρ - r2)) * (A + B + r1 - r2) / (4 * A * B)
    K = _carlson_rf(zero(m1), m1, one(m1))
    norm = sqrt(-_e2m1(E) * A * B)
    # from r1 down to r1 − d1 = r2 + d2: θ = π/2 + asin y, sin²θ = 1 − y², cos θ = −y, and
    # F(θ) = 2K − F(π − θ) above π/2
    function λ(d1, d2)
        den = B * d1 + A * d2
        y = (B * d1 - A * d2) / den
        f = _legendre_f(4 * A * B * d1 * d2 / den^2, y^2, m1)
        return (y <= 0 ? f : 2K - f) / norm
    end
    function λ_of_r(r)
        r2 <= r <= r1 && return λ(r1 - r, r - r2)
        @info("r = $r is out of the plunge region between r2 = $r2 and r1 = $r1.")
    end
    Λr_max = λ_of_r(r2)
    d1, d2 = gaps
    (d1 >= 0 && d2 >= 0) || return Λr_max, λ_of_r(_rplus(a)), λ_of_r   # r₊ outside the range
    λH = λ(d1, d2)
    return Λr_max, λH, _HorizonMinoMap(λ_of_r, _rplus(a), λH)
end

function lambda_of_r_real2(a, E, roots; gaps=(roots[4] - _rplus(a), _rplus(a) - roots[3]))
    r4, r3, r2, r1 = roots
    rp = _rplus(a)
    m1 = (r2 - r3) * (r1 - r4) / ((r1 - r3) * (r2 - r4))
    ξr = sqrt(-_e2m1(E) * (r1 - r3) * (r2 - r4)) / 2
    # from r1 down to r1 − d1 = r2 + d2: K − F(φ) = F(ψ), tan ψ = cot φ / k′, with
    # sin²φ = u = d2 (r1 − r3)/((r1 − r2)(r − r3)) and 1 − u = d1 (r2 − r3)/((r1 − r2)(r − r3))
    function λ(d1, d2)
        den = (r1 - r2) * (r2 + d2 - r3)
        u = d2 * (r1 - r3) / den; v = d1 * (r2 - r3) / den; w = v + u * m1
        return _legendre_f(v / w, u * m1 / w, m1) / ξr
    end
    function λ_of_r(r)
        rp <= r <= r1 && return λ(r1 - r, r - r2)
        @info("r = $r is out of the Real2 exterior plunge region between r+ = $rp and r1 = $r1.")
    end
    d1, d2 = gaps
    (d1 >= 0 && d2 >= 0) || return λ_of_r(rp), λ_of_r(rp), λ_of_r     # r₊ outside the range
    Λr_max = λ(d1, d2)
    return Λr_max, Λr_max, _HorizonMinoMap(λ_of_r, rp, Λr_max)
end

"""
    lambda_of_r(a, E, L, Q)

Return `(λ_end, λ_H, λ_of_r)` for the root class of `classify_orbit`: the Mino times from the
outer turning point to the inner end of the radial range (the inner turning point; r₊ for
Real2) and to the horizon r₊, and the map r ↦ λ measured from the outer turning point.
The outer turning point is r3 (Real1) or r1 (Real2, Complex). `λ_of_r(r)` is defined on the
radial range of the class; outside it logs the range with `@info` and returns `nothing`.
Root structures without an E < 1 plunge raise the error of `classify_orbit`.
"""
function lambda_of_r(a, E, L, Q)
    roots, cf = classify_orbit(a, E, L, Q)
    gaps(outer, inner) = (_horizon_gap(a, E, L, Q, outer), -_horizon_gap(a, E, L, Q, inner))
    if cf == "Real1"
        return lambda_of_r_real(a, E, roots; gaps=gaps(roots[2], roots[1]))
    elseif cf == "Complex"
        z = argmax(z -> abs(imag(z)), kerr_geo_root_structure(a, E, L, Q).raw_roots)
        return lambda_of_r_complex(a, E, roots; gaps=gaps(roots[1], roots[2]),
            pair=(real(z), abs(imag(z))))
    else   # "Real2", the only other class classify_orbit returns
        return lambda_of_r_real2(a, E, roots; gaps=gaps(roots[4], roots[3]))
    end
end
