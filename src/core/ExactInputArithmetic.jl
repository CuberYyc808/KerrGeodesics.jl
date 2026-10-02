# Error-free Float64 expansions certify zero only when every product residual is representable.
module _InputArithmetic

const Expansion = Vector{Float64}
struct Uncertified <: Exception end
scalar(x) = iszero(x) ? Float64[] : [Float64(x)]

function grow(e::Expansion, b::Float64)
    out = Float64[]
    total = b
    for a in e
        value = total + a
        isfinite(value) || throw(Uncertified())
        remainder = value - total
        low = (total - (value - remainder)) + (a - remainder)
        iszero(low) || push!(out, low)
        total = value
    end
    iszero(total) || push!(out, total)
    return out
end

function add(a::Expansion, b::Expansion)
    out = copy(a)
    for value in b
        out = grow(out, value)
    end
    return out
end
sub(a::Expansion, b::Expansion) = add(a, -b)

function mul(a::Expansion, b::Expansion)
    out = Float64[]
    for x in a, y in b
        # Bound the complete product of two 53-bit significands, including its low bits.
        exponent(x) + exponent(y) - 104 >= -1074 || throw(Uncertified())
        value = x * y
        isfinite(value) || throw(Uncertified())
        out = grow(grow(out, fma(x, y, -value)), value)
    end
    return out
end

function radial_coefficients(a, E, L, Q)
    aa, ee, ll, qq = scalar.((a, E, L, Q))
    a2 = mul(aa, aa)
    k = mul(sub(ee, scalar(1.0)), add(ee, scalar(1.0)))
    w = sub(mul(aa, ee), ll)
    return [-mul(a2, qq), mul(scalar(2.0), add(mul(w, w), qq)),
        sub(mul(a2, k), add(mul(ll, ll), qq)), scalar(2.0), k]
end

derivative(p) = length(p) == 1 ? [Float64[]] :
    [mul(scalar(k), p[k+1]) for k in 1:length(p)-1]

function evaluate(p, r)
    value = Float64[]
    rr = scalar(r)
    for coefficient in reverse(p)
        value = add(mul(value, rr), coefficient)
    end
    return value
end

function repeated_at(a, E, L, Q, r, multiplicity)
    try
        p = radial_coefficients(a, E, L, Q)
        for _ in 1:multiplicity
            isempty(evaluate(p, r)) || return false
            p = derivative(p)
        end
        return true
    catch err
        err isa Uncertified || rethrow()
        return false
    end
end

function polar_discriminant(a, E, L, Q)
    aa, ee, ll, qq = scalar.((a, E, L, Q))
    beta = mul(mul(aa, aa), mul(sub(ee, scalar(1.0)), add(ee, scalar(1.0))))
    b = sub(sub(beta, qq), mul(ll, ll))
    return add(mul(b, b), mul(scalar(4.0), mul(beta, qq)))
end

end
