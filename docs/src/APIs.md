# API Reference

This page documents public interfaces intended for direct user calls.

## Unified Constructors

```@docs
kerr_geodesic
```

```@docs
kerr_geo_stable
```

```@docs
kerr_geo_plunge
```

## Structured Output Types

```@docs
KerrGeodesicFamily
```

```@docs
KerrGeoStable
```

```@docs
KerrGeoPlunge
```

## Stable-Orbit API

Stable bound-orbit helper functions use the lower-case `kerr_geo_*` naming
scheme.

| function | purpose |
| :--- | :--- |
| `kerr_geo_constants_of_motion(a,p,e,x)` | return `Dict("E"=>..., "Lz"=>..., "Q"=>...)` |
| `kerr_geo_four_velocity(a,p,e,x; kwargs...)` | return stable-orbit four-velocity callables |
| `kerr_geo_frequencies(a,p,e,x; Time="Mino")` | return Mino, Boyer-Lindquist, or proper-time frequencies |
| `kerr_geo_orbit_type(a,p,e,x)` | return orbit-type labels |
| `kerr_geo_orbit(a,p,e,x; initPhases=...)` | return the dictionary-style stable geodesic object |

```@docs
kerr_geo_orbit_type_metadata
```

`kerr_geo_orbit_type_metadata` is the preferred source of stable-orbit
classification fields. `kerr_geo_orbit_type` returns labels derived from that
metadata; near separatrix inputs inside the roundoff guard are labeled
`Separatrix`.

## Plunge Root and Trajectory Helpers

```@docs
kerr_rstar
```

| function | purpose |
| :--- | :--- |
| `radial_roots(a,E,L,Q)` | return the radial-potential roots |
| `polar_roots(a,E,L,Q)` | return the polar-sector roots |
| `classify_orbit(a,E,L,Q; atol=1e-15)` | classify radial roots as `Complex`, `Real1`, `Real2`, or unsupported |
| `lambda_of_r(a,E,L,Q)` | return branch support radii and a radius-to-Mino-time map |
| `generic_plunge_orbit(a,E,L,Q; initPhases=...)` | return low-level plunge trajectory callables |
| `generic_plunge_velocity(a,E,L,Q; initPhase=...)` | return low-level plunge velocity callables |

The direct public plunge entry point is `kerr_geo_plunge`. Internal near-horizon time helpers are housed under
`KerrGeoPlunge/NearHorizonTime.jl`, and internal plunge initial-condition
helpers under `KerrGeoPlunge/InitialConditions.jl`; they are not separate public
constructor layers.

## Special-Orbit Helpers

| function | purpose |
| :--- | :--- |
| `kerr_geo_separatrix(a,e,x)` | separatrix radius for bound geodesics |
| `kerr_geo_isco(a,x)` | innermost stable circular orbit where implemented |
| `kerr_geo_ibso(a,x)` | innermost bound spherical orbit |
| `kerr_geo_isso(a,x)` | innermost stable spherical orbit |
