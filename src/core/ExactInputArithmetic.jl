# Error-free expansions (in the binary, correctly rounded floating-point type T of the
# constants, which has fma) certify zero only when every product residual is representable.
module _InputArithmetic

const Expansion{T} = Vector{T}
struct Uncertified <: Exception end
scalar(::Type{T}, x) where {T} = iszero(x) ? T[] : [T(x)]

function grow(e::Expansion{T}, b::T) where {T}
    out = T[]
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

function add(a::Expansion{T}, b::Expansion{T}) where {T}
    out = copy(a)
    for value in b
        out = grow(out, value)
    end
    return out
end
sub(a::Expansion{T}, b::Expansion{T}) where {T} = add(a, -b)

function mul(a::Expansion{T}, b::Expansion{T}) where {T}
    out = T[]
    p = precision(T)
    for x in a, y in b
        # Bound the complete product of two p-bit significands, including its low bits
        # (Float64: exponent sum − 104 ≥ −1074).
        exponent(x) + exponent(y) - 2(p - 1) >= exponent(floatmin(T)) - (p - 1) ||
            throw(Uncertified())
        value = x * y
        isfinite(value) || throw(Uncertified())
        out = grow(grow(out, fma(x, y, -value)), value)
    end
    return out
end

function radial_coefficients(a, E, L, Q)
    T = float(promote_type(typeof(a), typeof(E), typeof(L), typeof(Q)))
    aa, ee, ll, qq = scalar.(T, (a, E, L, Q))
    a2 = mul(aa, aa)
    k = mul(sub(ee, scalar(T, 1)), add(ee, scalar(T, 1)))
    w = sub(mul(aa, ee), ll)
    return [-mul(a2, qq), mul(scalar(T, 2), add(mul(w, w), qq)),
        sub(mul(a2, k), add(mul(ll, ll), qq)), scalar(T, 2), k]
end

derivative(p::Vector{Expansion{T}}) where {T} = length(p) == 1 ? [T[]] :
    [mul(scalar(T, k), p[k+1]) for k in 1:length(p)-1]

function evaluate(p::Vector{Expansion{T}}, r) where {T}
    value = T[]
    rr = scalar(T, r)
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
    T = float(promote_type(typeof(a), typeof(E), typeof(L), typeof(Q)))
    aa, ee, ll, qq = scalar.(T, (a, E, L, Q))
    beta = mul(mul(aa, aa), mul(sub(ee, scalar(T, 1)), add(ee, scalar(T, 1))))
    b = sub(sub(beta, qq), mul(ll, ll))
    return add(mul(b, b), mul(scalar(T, 4), mul(beta, qq)))
end

end
