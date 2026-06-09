module OrbitClass

using Roots
using Polynomials

export radial_roots, polar_roots, classify_orbit

"""
    radial_roots(a, E, L, Q)

Return the four radial-potential roots for a Kerr geodesic with constants of
motion `(E,L,Q)`.
"""
function radial_roots(a, E, L, Q)
    # Solve the fourth-order polynomial
    coe4 = E^2 - 1
    coe3 = 2
    coe2 = - (Q + L^2 + a^2 * (1 - E^2))
    coe1 = 2 * (a * E - L)^2 + 2 * Q
    coe0 = - a^2 * Q
    radial_zeros = roots(Polynomial([coe0, coe1, coe2, coe3, coe4]))
    return radial_zeros
end

"""
    polar_roots(a, E, L, Q)

Return the polar-sector root pair used by the plunge trajectory layer.
"""
function polar_roots(a, E, L, Q)
    A = (Q + L^2) / (a^2 * (1 - E^2)) + 1
    B = Q / (a^2 * (1 - E^2))
    zm = (A - sqrt(A^2 - 4 * B)) / 2
    zp = (A + sqrt(A^2 - 4 * B)) / 2
    return zm, zp
end

"""
    classify_orbit(a, E, L, Q; atol=1e-15)

Classify the radial-root structure as `Complex`, `Real1`, `Real2`, or
`Unsupported`.
"""
function classify_orbit(a, E, L, Q; atol=1e-15)
    # Compute radial roots
    roots = radial_roots(a, E, L, Q)
    rp = 1 + sqrt(1 - a^2)

    real_roots = Float64[]
    complex_roots = ComplexF64[]

    for r in roots
        if abs(imag(r)) < atol
            push!(real_roots, real(r))
        else
            push!(complex_roots, r)
        end
    end

    if length(real_roots) == 4
        sort!(real_roots)

        n_outside = count(r -> r > rp, real_roots)

        if n_outside == 3
            return real_roots, "Real1"
        elseif n_outside == 1
            return real_roots, "Real2"
        else
            error("Unexpected real-root configuration relative to r_+")
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
            "Unexpected root structure: " *
            "$(length(real_roots)) real roots and $(length(complex_roots)) complex roots"
        )
    end
end

end
