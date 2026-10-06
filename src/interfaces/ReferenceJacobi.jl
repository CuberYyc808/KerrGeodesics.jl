# Jacobi amplitude and sn, cn, dn for the APEX and plunge reference interfaces, whose numbers are kept unchanged.
# The arithmetic is that of Elliptic.jl's Jacobi module (MIT license): the Landen sequence of Abramowitz & Stegun
# 16.4 for 0 <= m <= 1, the series 16.13.4 and 16.15.4 next to m = 0 and m = 1, and 16.10 for m outside [0, 1].
# The values are therefore bitwise the same. The only difference is where the Landen ratios c_n/a_n are held:
# Elliptic.jl keeps them in a module-level buffer that concurrent calls overwrite, here they are a tuple local
# to each call, so trajectories can be evaluated from several threads at once.

function _reference_am(u::Float64, m::Float64)
    u == 0.0 && return 0.0
    tol = eps(Float64)
    sqrt_tol = sqrt(tol)
    m < sqrt_tol && return u - 0.25 * m * (u - 0.5 * sin(2.0 * u))
    m1 = 1.0 - m
    if m1 < sqrt_tol
        t = tanh(u)
        return asin(t) + 0.25 * m1 * (t - u * (1.0 - t^2)) * cosh(u)
    end
    ratios = ntuple(_ -> 0.0, Val(10))
    a, b, c, n = 1.0, sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        n < 10 || error("Landen sequence did not converge in 10 steps for m = $m")
        a, b, c, n = 0.5 * (a + b), sqrt(a * b), 0.5 * (a - b), n + 1
        ratios = Base.setindex(ratios, c / a, n)
    end
    phi = ldexp(a * u, n)
    for i in n:-1:1
        phi = 0.5 * (phi + asin(ratios[i] * sin(phi)))
    end
    return phi
end

function _reference_am_checked(u::Real, m::Real)
    (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
    return _reference_am(Float64(u), Float64(m))
end

# sn, cn, dn: Abramowitz & Stegun 16.10 maps m < 0 and m > 1 onto the amplitude at a parameter in [0, 1]
function _reference_sncndn(u::Float64, m::Float64)
    if m < 0.0
        mu1 = 1.0 / (1.0 - m)
        mu = -m * mu1
        sqrtmu1 = sqrt(mu1)
        phi = _reference_am(u / sqrtmu1, mu)
        s = sin(phi)
        d = sqrt(1.0 - mu * s^2)
        return (sqrtmu1 * s) / d, cos(phi) / d, 1.0 / d
    elseif m > 1.0
        mu = 1 / m
        phi = _reference_am(u * sqrt(m), mu)
        return sqrt(mu) * sin(phi), sqrt(1.0 - mu * sin(phi)^2), cos(phi)
    end
    phi = _reference_am(u, m)
    return sin(phi), cos(phi), sqrt(1.0 - m * sin(phi)^2)
end

_reference_sn(u::Real, m::Real) = _reference_sncndn(Float64(u), Float64(m))[1]
_reference_cn(u::Real, m::Real) = _reference_sncndn(Float64(u), Float64(m))[2]
_reference_dn(u::Real, m::Real) = _reference_sncndn(Float64(u), Float64(m))[3]
