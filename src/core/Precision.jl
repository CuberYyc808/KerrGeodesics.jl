# Working precision. Every computation runs in the floating-point type T of its inputs
# (`_float_type`); tolerances and truncation orders chosen for Float64 are carried to T by
# `_tol` and `_nterms`, which return the Float64 values unchanged. An orbit built in BigFloat
# remembers the precision it was built with and evaluates under it (`_with_precision`).

# the floating-point type of a computation on these inputs
_float_type(xs...) = float(promote_type(map(typeof, xs)...))
# the real floating-point type of real or complex inputs
_real_type(xs...) = real(_float_type(xs...))

# An empirical Float64 tolerance sits between rounding (eps) and the physical scale (1): it
# keeps its place on the logarithmic scale between them, tol(T) = tol^(log eps(T)/log eps(Float64)).
_tol(::Type{T}, tol) where {T} = T(tol)^(log(eps(T)) / log(eps(Float64)))
_tol(::Type{Float64}, tol) = tol

# the number of terms of a series truncated after n64 terms in Float64, for the precision of T
# (a geometric series needs a number of terms proportional to the number of digits)
_nterms(::Type{T}, n64) where {T} = cld(n64 * precision(T), 53)

# evaluate f() with BigFloat at precision p (the precision an object was built with); other
# types carry their precision in the type itself
_with_precision(f, ::Type{BigFloat}, p) = setprecision(f, BigFloat, p)
_with_precision(f, ::Type, p) = f()

# the precision of BigFloat inputs (their largest), under which a result is built
_input_precision(xs...) = maximum(x -> x isa BigFloat ? precision(x) : 0, xs; init=0)

# A function of a BigFloat result evaluated, wherever it is called, at the precision p the
# result was built with; other types need no wrapping
_precision_wrap(::Type{BigFloat}, p, f::Function) =
    (args...; kwargs...) -> setprecision(() -> f(args...; kwargs...), BigFloat, p)
_precision_wrap(::Type{BigFloat}, p, x::Union{NamedTuple,Tuple}) =
    map(v -> _precision_wrap(BigFloat, p, v), x)
_precision_wrap(::Type{BigFloat}, p, x::AbstractVector) = [_precision_wrap(BigFloat, p, v) for v in x]
_precision_wrap(::Type{BigFloat}, p, x) = x
_precision_wrap(::Type, p, x) = x
