# Every floating-point number reachable from a result (through tuples, named tuples, arrays,
# dictionaries and struct fields), with its path: a BigFloat computation must return no Float64.
function _floats!(out, x, path, seen)
    if x isa AbstractFloat
        push!(out, (path, typeof(x)))
    elseif x isa Complex
        _floats!(out, real(x), path * ".re", seen); _floats!(out, imag(x), path * ".im", seen)
    elseif x isa Union{Tuple,NamedTuple}
        for (k, v) in pairs(x)
            _floats!(out, v, "$path.$k", seen)
        end
    elseif x isa AbstractArray
        for (k, v) in enumerate(x)
            _floats!(out, v, "$path[$k]", seen)
        end
    elseif x isa AbstractDict
        for (k, v) in x
            _floats!(out, v, "$path[$k]", seen)
        end
    elseif x isa Union{Function,Symbol,AbstractString,Number,Nothing,Type,Module,Bool}
    elseif isstructtype(typeof(x)) && !(x in seen)
        push!(seen, x)
        for k in fieldnames(typeof(x))
            isdefined(x, k) && _floats!(out, getfield(x, k), "$path.$k", seen)
        end
    end
    return out
end
float_leaves(x) = _floats!(Tuple{String,DataType}[], x, "", IdSet{Any}())
non_bigfloat(x) = [p for (p, t) in float_leaves(x) if t !== BigFloat]
