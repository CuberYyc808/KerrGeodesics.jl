# Catalogue constants in a floating-point type T. Most catalogue cases are generic (binary64
# constants, exact in every precision). A repeated root, a root on the horizon or a triple
# polar-coefficient condition is a property of the exact constants, so for those cases the
# remaining constants are solved in T: P(r₊) = 0 by Lz = 2r₊E/a, h = 0 at |a| = 1 by
# Q = 3E² − 1 (or Q = E⁴/(1 − E²) for A-X2), and a radial root of multiplicity m by Newton's
# method on R = R′ = … = R^(m−1) = 0 for the root and m − 1 constants, from the binary64 values.

# R(r) = Σ c_j r^j and its k-th r-derivative at r
_rcoeffs(a, E, L, Q) = (-a^2 * Q, 2 * ((a * E - L)^2 + Q), a^2 * (E^2 - 1) - L^2 - Q, 2one(E), E^2 - 1)
function _Rk(k, r, a, E, L, Q)
    c = _rcoeffs(a, E, L, Q)
    return sum(c[j + 1] * (factorial(j) ÷ factorial(j - k)) * r^(j - k) for j in k:4)
end

# Newton on R = R′ = … = R^(m−1) = 0 for x = (r, unknowns...), with g(x) = (a, E, L, Q);
# partial derivatives by complex steps (no cancellation)
function _solve_repeated(::Type{T}, g, x0, m) where {T}
    x = collect(T, x0)
    h = T(2)^(-2precision(T))
    F(y) = [_Rk(k, y[1], g(y)...) for k in 0:m-1]
    for _ in 1:100
        J = zeros(T, m, m)
        for j in 1:m
            y = complex.(x); y[j] += im * h
            J[:, j] = imag.(F(y)) ./ h
        end
        step = J \ F(x)
        x -= step
        maximum(abs.(step) ./ max.(abs.(x), 1)) <= 4eps(T) && break
    end
    return x
end

function exact_catalogue_constants(::Type{T}, id, a, E, L, Q) where {T}
    a, E, L, Q = T(a), T(E), T(L), T(Q)
    rp = 1 + sqrt(1 - a^2)
    id in (:A_H1, :D_H1, :D_H2) && return (a, E, 2rp * E / a, Q)
    # A-X2: Q = E⁴/(1 − E²) rounded down, so that the extremal quadratic's discriminant is not
    # rounded below zero (the stable circular orbit is approached from the side of real roots)
    if id === :A_X2
        T === Float64 && return (a, E, L, E^4 / (1 - E^2))
        q = setprecision(BigFloat, 2precision(T)) do
            e = BigFloat(E); e^4 / (1 - e^2)
        end
        return (a, E, L, BigFloat(q, RoundDown; precision=precision(T)))
    end
    id in (:B_X2, :C_X4) && return (a, E, L, 3E^2 - 1)
    rootdata = Dict(:A2 => (8.0, 2), :B2 => (8.0, 2), :K3 => (3.2265812774414178, 2),
        :K4 => (3.2265812774414178, 2), :K5 => (3.2265812774414178, 2), :K6 => (4.0, 2),
        :K7 => (4.0, 2), :K8 => (4.0, 2), :K9 => (3.0, 2), :K10 => (3.0, 2), :K11 => (3.0, 2),
        :B7 => (0.01, 2), :B8 => (0.05, 2), :C6 => (0.01, 2), :C7 => (0.05, 2),
        :C9 => (0.01, 2), :C10 => (0.07, 2), :C11 => (-0.28, 2), :N2 => (8.0, 2),
        :K1 => (2.794219261817521, 3), :K2 => (2.794219261817521, 3),
        :A_H2 => (7.015979675282077, 2))
    haskey(rootdata, id) || return (a, E, L, Q)
    r0, m = rootdata[id]
    if id === :A_H2           # Q = 0, Lz = 2r₊E/a: solve for (r, E)
        x = _solve_repeated(T, y -> (a, y[2], 2rp * y[2] / a, Q), (r0, E), 2)
        return (a, x[2], 2rp * x[2] / a, Q)
    elseif m == 3             # triple root: solve for (r, Lz, Q) at fixed E
        x = _solve_repeated(T, y -> (a, E, y[2], y[3]), (r0, L, Q), 3)
        return (a, E, x[2], x[3])
    end
    # (r, Lz) at fixed E, Q; where ∂R/∂Lz vanishes at the root (K6–K8), (r, Q) at fixed E, Lz
    if id in (:K6, :K7, :K8)
        x = _solve_repeated(T, y -> (a, E, L, y[2]), (r0, Q), 2)
        return (a, E, L, x[2])
    end
    x = _solve_repeated(T, y -> (a, E, y[2], Q), (r0, L), 2)
    return (a, E, x[2], Q)
end
