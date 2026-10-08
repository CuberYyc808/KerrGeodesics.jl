# Plunge reference API: radial and polar roots and the Real1/Real2/Complex root classes of an E < 1 plunge.

"""
    radial_roots(a, E, L, Q)

The four roots of the radial potential R(r) of the constants `(E, L, Q)` (L = Lz), as a
vector of complex numbers.
"""
function radial_roots(a, E, L, Q)
    # Solve the fourth-order polynomial
    coe4 = E^2 - 1
    coe3 = 2
    coe2 = - (Q + L^2 + a^2 * (1 - E^2))
    coe1 = 2 * (a * E - L)^2 + 2 * Q
    coe0 = - a^2 * Q
    radial_zeros = _polynomial_roots([coe0, coe1, coe2, coe3, coe4])
    return radial_zeros
end

"""
    polar_roots(a, E, L, Q)

Return `(zm, zp)`, the roots zm ≤ zp of a²(1 − E²) y² − (Q + L² + a²(1 − E²)) y + Q = 0 in
y = z² = cos²θ; zp = Inf when a²(1 − E²) = 0.
"""
function polar_roots(a, E, L, Q)
    c = a^2 * (1 - E^2)
    roots = _polar_quadratic_roots(c, L, Q)
    return roots.u_small, c > 0 ? roots.u_big : oftype(roots.u_small, Inf)
end

"""
    _plunge_polar_parameters(a, E, L, Q) -> (zm, ξθ, kθ)

`z = sqrt(zm) sn(ξθ (λ + λθ0), kθ)` with `ξθ^2 = a^2(1-E^2) zp` and `kθ = zm/zp`, written so
that the a -> 0 limit (ξθ^2 = Q + L^2, kθ = 0) is exact.
"""
function _plunge_polar_parameters(a, E, L, Q)
    c = a^2 * (1 - E^2)
    roots = _polar_quadratic_roots(c, L, Q)
    return roots.u_small, sqrt(roots.cu_big), roots.cu_small / roots.cu_big
end

"""
    classify_orbit(a, E, L, Q; atol=1e-15)

Return `(roots, class)` for an E < 1 plunge: `"Real1"` (four real roots, three outside r₊) or
`"Real2"` (four real roots, one outside r₊) with `roots` ascending, r4 ≤ r3 ≤ r2 ≤ r1;
`"Complex"` (two real roots r2 < r1 and a complex pair ρ, ρ̄) with `roots = [r1, r2, A, B]`,
A = |r1 − ρ|, B = |r2 − ρ|. Any other root structure is an error.
The keyword `atol` is retained for call compatibility; root reality is determined by
conjugate pairing of the refined roots.
"""
function classify_orbit(a, E, L, Q; atol=1e-15)
    iszero(_wide_horizon_momentum(a, E, L)) && error(
        "classify_orbit: a plunge with zero horizon momentum requires a horizon-root formula.")
    # the four roots refined together in double-double (`kerr_geo_root_structure`); these
    # closed forms need four distinct roots, so no repeated-root reading is applied. A root is
    # real when the estimate nearest its conjugate is itself, nonreal when it is another one.
    raw = collect(kerr_geo_root_structure(a, E, L, Q).raw_roots)
    rp = _rplus(a)
    paired = [argmin(w -> abs(w - conj(z)), raw) != z for z in raw]
    # nonreal roots of a real polynomial come in pairs: with an odd count (a cluster of three
    # nearly equal roots), the one nearest the axis is real
    if isodd(count(paired))
        paired[argmin(i -> paired[i] ? abs(imag(raw[i])) : Inf, eachindex(raw))] = false
    end
    T = real(eltype(raw))
    real_roots = T[real(raw[i]) for i in eachindex(raw) if !paired[i]]
    complex_roots = Complex{T}[raw[i] for i in eachindex(raw) if paired[i]]

    if length(real_roots) == 4
        sort!(real_roots)

        # R(r₊) = P(r₊)² > 0 puts r₊ inside an allowed interval, so an odd number of roots lie
        # outside it; an even count means the root nearest r₊ rounded across it (P(r₊) → 0)
        n_outside = count(r -> r > rp, real_roots)
        if iseven(n_outside) && !iszero(_wide_horizon_momentum(a, E, L))
            n_outside += argmin(r -> abs(r - rp), real_roots) > rp ? -1 : 1
        end

        if n_outside == 3
            return real_roots, "Real1"
        elseif n_outside == 1
            return real_roots, "Real2"
        else
            error("classify_orbit: $n_outside of the four real roots lie outside r₊; an E < 1 " *
                "plunge has three (Real1) or one (Real2).")
        end

    elseif length(real_roots) == 2 && length(complex_roots) == 2
        sort!(real_roots)

        r2 = real_roots[1]   # smaller
        r1 = real_roots[2]   # larger

        ρr = real(complex_roots[1])
        ρi = abs(imag(complex_roots[1]))

        # Define auxiliary quantities
        A = sqrt((r1 - ρr)^2 + ρi^2)
        B = sqrt((r2 - ρr)^2 + ρi^2)

        return [r1, r2, A, B], "Complex"
    else
        error(
            "classify_orbit: $(length(real_roots)) real and $(length(complex_roots)) complex " *
            "roots; an E < 1 plunge has four real roots (Real1, Real2) or two real roots and a " *
            "complex pair (Complex)."
        )
    end
end
