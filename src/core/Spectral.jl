# Adaptive piecewise Chebyshev representation of smooth functions, and of their integrals.
#
# The Boyer–Lindquist t and φ and the proper time τ of a Kerr geodesic are Mino-time
# integrals of rates that split into a part depending on r alone and a part depending on z
# alone. The rates are sampled once along the analytic
# r(λ), z(λ) and stored as Chebyshev series on adaptively chosen pieces; integrating a
# Chebyshev series is exact, so the primitive is available to the same relative accuracy
# and evaluates in tens of nanoseconds (find the piece, one Clenshaw sum).

"""
Piecewise Chebyshev series of `ncomp` functions on the pieces `breaks[i]..breaks[i+1]`.
`achieved[i]` is the error estimate of piece i relative to the local size of the function
(worst component): the coefficient tail or the misfit at off-node points, whichever is
larger.
"""
struct ChebPieces
    breaks::Vector{Float64}
    coefs::Vector{Vector{Vector{Float64}}}      # coefs[component][piece]
    achieved::Vector{Float64}
end

const _CHEB_N = 32

_cheb_nodes(n) = n == _CHEB_N ? _CHEB_NODES : [cospi(j / n) for j in 0:n]  # x_j on [-1, 1], x_0 = 1
const _CHEB_NODES = [cospi(j / _CHEB_N) for j in 0:_CHEB_N]

# cos(π j k / n), tabulated for the default order (the DCT below is the build cost of
# every table: n² cosines per component and piece otherwise)
_dct_table(n) = [cospi(j * k / n) for j in 0:n, k in 0:n]
const _DCT_TABLE = _dct_table(_CHEB_N)

# DCT-I: values at the n+1 Chebyshev points -> coefficients of Σ c_k T_k
function _cheb_coefficients(values::AbstractVector{Float64})
    n = length(values) - 1
    C = n == _CHEB_N ? _DCT_TABLE : _dct_table(n)
    c = zeros(n + 1)
    @inbounds for k in 0:n
        s = 0.5 * (values[1] + (isodd(k) ? -1 : 1) * values[n + 1])
        for j in 1:n-1
            s += values[j + 1] * C[j + 1, k + 1]
        end
        c[k + 1] = 2s / n
    end
    c[1] /= 2
    c[n + 1] /= 2
    return c
end

@inline function _clenshaw(c::Vector{Float64}, x::Float64)
    b1 = 0.0; b2 = 0.0
    @inbounds for k in length(c):-1:2
        b1, b2 = muladd(2x, b1, c[k] - b2), b1
    end
    return muladd(x, b1, c[1] - b2)
end

# Divided Clenshaw recurrence: (f(x)-f(y))/(x-y), without subtracting
# two primitives. The physical x-y is supplied separately after rescaling.
@inline function _clenshaw_divided(c::Vector{Float64}, x::Float64, y::Float64)
    b1=0.0; b2=0.0; d1=0.0; d2=0.0
    @inbounds for k in length(c):-1:2
        d1,d2=muladd(2x,d1,2b1-d2),d1
        b1,b2=muladd(2y,b1,c[k]-b2),b1
    end
    return muladd(x,d1,b1-d2)
end

# Off-node check points (in [-1, 1], away from the n = 32 Chebyshev points): a function that
# the nodes alias (e.g. T_64 on 33 points looks constant) or a feature narrower than the
# node spacing shows up as a misfit there.
const _CHEB_CHECK = (-0.8763, -0.3371, 0.2197, 0.6639)

# One piece: sample, and return the chopped coefficients, the excess and the achieved
# accuracy. The excess is the largest ratio (error estimate) / (allowed error) over the
# components, where the error estimate is the larger of the coefficient tail and the
# off-node misfit, and the allowed error is tol × local size (or `absfloor[k]`), or the
# abscissa-rounding noise, or the rounding of the values themselves `abserr[k]`, whichever is
# larger. The piece is resolved when the excess is ≤ 1. `nothing` for non-finite samples.
function _fit_piece(f, a, b, ncomp, tol, absfloor, abserr, n)
    xs = _cheb_nodes(n)
    vals = Matrix{Float64}(undef, n + 1, ncomp)
    for (j, x) in enumerate(xs)
        # (the end nodes are the interval ends exactly: f may be defined on [a, b] only)
        v = f(j == 1 ? b : j == n + 1 ? a : (a + b) / 2 + (b - a) / 2 * x)
        for k in 1:ncomp
            vk = v[k]
            isfinite(vk) || return nothing
            vals[j, k] = vk
        end
    end
    checks = map(x -> f((a + b) / 2 + (b - a) / 2 * x), _CHEB_CHECK)
    # the value-rounding floor of this piece: a constant, or a function of x (the larger of its
    # values at the two ends; it varies slowly along a piece)
    floorv = abserr isa Function ? max.(abserr(a), abserr(b)) : abserr
    out = Vector{Vector{Float64}}(undef, ncomp)
    excess = 0.0
    achieved = 0.0
    for k in 1:ncomp
        c = _cheb_coefficients(view(vals, :, k))
        own = maximum(abs, view(vals, :, k))
        scale = max(own, absfloor[k])
        # Rounding of the abscissa itself limits what any piece can resolve: a relative
        # error eps in x moves f by eps·|x|·|f'|. Where f is that steep the tail is noise.
        slope = 0.0
        for j in 2:n+1
            slope += (j - 1)^2 * abs(c[j])
        end
        noise = 8eps() * max(abs(a), abs(b)) * slope * 2 / (b - a)
        bound = max(tol * scale, noise, floorv[k], floatmin())
        tail = maximum(abs, view(c, n-2:n+1))
        misfit = 0.0
        for (x, v) in zip(_CHEB_CHECK, checks)
            isfinite(v[k]) || return nothing
            misfit = max(misfit, abs(_clenshaw(c, x) - v[k]))
        end
        estimate = max(tail, misfit)
        # noise plateau: the coefficients stopped decaying far below the component's own
        # size (rounding in f itself; bisection cannot improve on it). Measured against the
        # component's own values, not the floor: a small but unresolved feature must not
        # pass for rounding noise.
        plateau = estimate <= 1.0e-9 * own && maximum(abs, view(c, n-14:n-7)) <= 4tail
        excess = max(excess, plateau ? min(1.0, estimate / bound) : estimate / bound)
        achieved = max(achieved, estimate / max(scale, floatmin()))
        cut = max(bound, tail) / 8
        m = n + 1
        while m > 1 && abs(c[m]) <= cut
            m -= 1
        end
        out[k] = c[1:m]
    end
    return out, excess, achieved
end

"""
    chebfit(f, breaks; ncomp=1, tol=1e-14, absfloor, abserr, maxdepth=60, maxpieces=2000)

Fit `f(x)` (returning `ncomp` values) on each interval of `breaks`, bisecting every
interval until its error estimate (Chebyshev tail, and the misfit at four off-node points)
is below `tol` times the local size of the function (or `absfloor[k]`), or below `abserr`,
the rounding of the sampled values themselves where the caller knows it (a tuple, or a
function of x evaluated at a piece's two ends; the radial engine next to a horizon: the rates
are evaluated from a rounded radius). Non-finite samples
force a split. A piece whose coefficients have levelled off at a noise plateau far below
the function's own size (rounding in the function being sampled), or whose tail is at the
rounding noise of the abscissa, is kept as it is: bisection cannot improve it. A piece
still unresolved when the depth, the piece budget or the minimum width is exhausted raises
an error. `achieved` of the result records each piece's estimate.
"""
function chebfit(f, breaks::AbstractVector{<:Real}; ncomp::Int=1, tol=1.0e-14,
        absfloor=zeros(ncomp), abserr=zeros(ncomp), maxdepth::Int=60, maxpieces::Int=2000,
        n::Int=_CHEB_N)
    outb = Float64[float(breaks[1])]
    outc = [Vector{Float64}[] for _ in 1:ncomp]
    achieved = Float64[]
    stack = Tuple{Float64,Float64,Int}[]
    for i in length(breaks)-1:-1:1
        push!(stack, (float(breaks[i]), float(breaks[i + 1]), 0))
    end
    while !isempty(stack)
        a, b, depth = pop!(stack)
        fit = _fit_piece(f, a, b, ncomp, tol, absfloor, abserr, n)
        excess = fit === nothing ? Inf : fit[2]
        if excess > 1
            splittable = depth < maxdepth && length(outb) + length(stack) < maxpieces &&
                b - a > 64eps(max(abs(a), abs(b), 1.0))
            if splittable
                m = (a + b) / 2
                push!(stack, (m, b, depth + 1))
                push!(stack, (a, m, depth + 1))
                continue
            end
            fit === nothing && error("chebfit: non-finite samples on [$a, $b].")
            error("chebfit: unresolved on [$a, $b] (error estimate $(fit[3]) of the " *
                "function's size, tolerance $tol; depth $depth, $(length(outb) - 1) pieces).")
        end
        push!(outb, b)
        for k in 1:ncomp
            push!(outc[k], fit[1][k])
        end
        push!(achieved, fit[3])
    end
    return ChebPieces(outb, outc, achieved)
end
chebfit(f, a::Real, b::Real; kwargs...) = chebfit(f, [a, b]; kwargs...)

npieces(p::ChebPieces) = length(p.breaks) - 1

@inline function _piece_index(p::ChebPieces, x)
    br = p.breaks
    i = searchsortedlast(br, x)
    return clamp(i, 1, length(br) - 1)
end

"""Evaluate component `k` at `x` (clamped to the covered interval's end pieces)."""
@inline function (p::ChebPieces)(x::Real, k::Int=1)
    i = _piece_index(p, x)
    a = p.breaks[i]; b = p.breaks[i + 1]
    return _clenshaw(p.coefs[k][i], (2x - a - b) / (b - a))
end

"""Evaluate the first three components at `x` (type-stable fast path)."""
@inline _eval3(p::ChebPieces, x) = (i = _piece_index(p, x); a = p.breaks[i]; b = p.breaks[i + 1];
    t = (2x - a - b) / (b - a);
    (_clenshaw(p.coefs[1][i], t), _clenshaw(p.coefs[2][i], t), _clenshaw(p.coefs[3][i], t)))

function _cheb_increment(p::ChebPieces, left, right, k)
    left==right && return 0.0
    right<left && return -_cheb_increment(p,right,left,k)
    first=_piece_index(p,left); last=_piece_index(p,right)
    value=0.0
    for i in first:last
        a=p.breaks[i]; b=p.breaks[i+1]
        lo=i==first ? left : a; hi=i==last ? right : b
        x=(2hi-a-b)/(b-a); y=(2lo-a-b)/(b-a)
        value+=2*((hi-lo)/(b-a))*_clenshaw_divided(p.coefs[k][i],x,y)
    end
    return value
end

function _cheb_increment_delta(p::ChebPieces, left, delta, k)
    iszero(delta) && return 0.0
    right = left + delta
    lo, hi = minmax(left, right)
    first = _piece_index(p, lo); last = _piece_index(p, hi)
    value = 0.0
    for i in first:last
        a = p.breaks[i]; b = p.breaks[i + 1]
        if first == last
            x = (2left - a - b) / (b - a)
            dx = 2delta / (b - a)
            value += dx * _clenshaw_divided(p.coefs[k][i], x + dx, x)
        else
            xl = i == first ? lo : a; xr = i == last ? hi : b
            width = if delta > 0
                i == first ? b - left : i == last ? delta - (a - left) : b - a
            else
                i == first ? delta + (left - b) : i == last ? a - left : a - b
            end
            x = (2xr - a - b) / (b - a); y = (2xl - a - b) / (b - a)
            value += 2width / (b - a) * _clenshaw_divided(p.coefs[k][i], x, y)
        end
    end
    return value
end

# ∫ Σ c_k T_k dx on [-1, 1] as a series vanishing at x = -1
function _cheb_integral(c::Vector{Float64}, halfwidth)
    n = length(c)
    C = zeros(n + 1)
    ext(k) = k < n ? c[k + 1] : 0.0
    for k in 1:n
        C[k + 1] = (ext(k - 1) * (k == 1 ? 2 : 1) - ext(k + 1)) / (2k)
    end
    # value at x = -1 is Σ C_k (-1)^k; shift C_0 so that it vanishes
    C[1] = -sum(C[k + 1] * (isodd(k) ? -1 : 1) for k in 1:n)
    return C .* halfwidth
end

"""
    chebintegrate(p) -> ChebPieces

Primitives of every component, continuous across pieces and zero at `p.breaks[1]`.
"""
function chebintegrate(p::ChebPieces)
    ncomp = length(p.coefs)
    out = [Vector{Vector{Float64}}(undef, npieces(p)) for _ in 1:ncomp]
    for k in 1:ncomp
        offset = 0.0
        for i in 1:npieces(p)
            h = (p.breaks[i + 1] - p.breaks[i]) / 2
            C = _cheb_integral(p.coefs[k][i], h)
            C[1] += offset
            offset = sum(C)                              # value at x = +1
            out[k][i] = C
        end
    end
    return ChebPieces(copy(p.breaks), out, copy(p.achieved))
end

"""Worst achieved accuracy and piece count of several tables: `(achieved, pieces)`."""
_spectral_summary(tables) = (achieved=maximum((maximum(t.achieved; init=0.0) for t in tables); init=0.0),
    pieces=sum((npieces(t) for t in tables); init=0))

"""Total integral of component `k` over the covered interval."""
function chebtotal(prim::ChebPieces, k::Int=1)
    return sum(prim.coefs[k][end])
end

"""
Spectral accuracy record of a member (`Status.spectral`): `achieved` is the worst error
estimate of all its Chebyshev tables relative to the local size of the tabulated rate,
`pieces` their number, `unresolved` the pieces left unresolved (0: `chebfit` raises an
error instead of returning one; pieces limited by rounding noise count as resolved and show
in `achieved`). Tables built lazily (the radial engine, on the first t or φ call)
are built when the record is read.
"""
struct SpectralStatus{F}
    summaries::F                        # () -> iterable of (achieved, pieces)
end
function Base.getproperty(s::SpectralStatus, k::Symbol)
    k === :summaries && return getfield(s, :summaries)
    parts = getfield(s, :summaries)()
    k === :achieved && return maximum((p.achieved for p in parts); init=0.0)
    k === :pieces && return sum((p.pieces for p in parts); init=0)
    k === :unresolved && return 0
    throw(ArgumentError("SpectralStatus has the properties achieved, pieces, unresolved"))
end
Base.propertynames(::SpectralStatus) = (:achieved, :unresolved, :pieces)
Base.show(io::IO, s::SpectralStatus) =
    print(io, "(achieved = ", s.achieved, ", unresolved = 0, pieces = ", s.pieces, ")")

const _NO_TABLES = SpectralStatus(() -> ())
