# Public access to a member's coordinate engine and to the horizon-relative radial state of the
# radius models that cross the future horizon. The functions only forward to the internal
# methods (`_radial_state`, `_relative_rstar`), so their values are those the members use.

"""
    engine_coordinates(member) -> EngineCoordinates or nothing

The [`EngineCoordinates`](@ref) behind the trajectory closures of `member` (a component such
as [`KerrGeoPlungeComponent`](@ref)), or `nothing` for a member whose coordinates are not
given by a coordinate engine.

```julia
m = kerr_geo_plunge_component(0.9, 0.95, 2.0, 3.0)
c = engine_coordinates(m)
state = horizon_relative_state(c, 1.0)
```
"""
function engine_coordinates(member)
    hasproperty(member, :Trajectory) || return nothing
    return _find_engine_coordinates(getproperty(member, :Trajectory), 0)
end

# depth-limited search through the fields of the trajectory NamedTuple and the closures in it
function _find_engine_coordinates(x, depth)
    x isa EngineCoordinates && return x
    depth >= 7 && return nothing
    (x isa Function || x isa NamedTuple || x isa Tuple) || return nothing
    for i in 1:fieldcount(typeof(x))
        found = _find_engine_coordinates(getfield(x, i), depth + 1)
        found === nothing || return found
    end
    return nothing
end

"""
    horizon_relative_state(r_of, λ) -> (gap, velocity, chart) or nothing
    horizon_relative_state(c::EngineCoordinates, λ)

The radial state at Mino time `λ` of a radius function that crosses the future horizon,
measured from the horizon: `gap = r - r₊`, formed without subtracting `r₊` from `r` (so it
keeps its relative accuracy down to the horizon), `velocity = dr/dλ`, and `chart`, the radius
model, whose fields `horizon` (``r_+`` as an unevaluated sum `(hi, lo)`) and `separation`
(``r_+ - r_-``) give the radius `r = chart.horizon[1] + (chart.horizon[2] + gap)`.
[`horizon_relative_rstar`](@ref) gives the tortoise coordinate from the state.

Returns `nothing` when `r_of` is not horizon-relative (members that do not reach the future
horizon); evaluate `r_of(λ)` then. With an [`EngineCoordinates`](@ref) the radius function is
its `r_of` (see [`engine_coordinates`](@ref)).
"""
horizon_relative_state(r_of, λ) = _radial_state(r_of, λ)
horizon_relative_state(c::EngineCoordinates, λ) = _radial_state(c.r_of, λ)

"""
    horizon_relative_rstar(state)

The tortoise coordinate ``r_*`` of a [`horizon_relative_state`](@ref), with the convention of
[`kerr_rstar`](@ref), computed from `state.gap`: the logarithm of ``r - r_+`` is taken of the
gap itself, so ``r_*`` stays accurate as the orbit approaches the horizon, where
`kerr_rstar(a, r)` loses the digits of ``r - r_+``. At ``r_+ = r_-`` it is the extremal form
``r + 2\\ln(r - 1) - 2/(r - 1) - 2\\ln 2``.
"""
horizon_relative_rstar(state) = _relative_rstar(state)
