# KerrGeodesics.jl

KerrGeodesics.jl provides Julia interfaces for Kerr geodesic trajectories in
units with $G=c=M=1$.

The package contains two project-facing layers:

- stable bound geodesics from APEX-like parameters `(a,p,e,x)`;
- bound plunge geodesics from constants of motion `(a,E,Lz,Q)`.

The unified constructor is [`kerr_geodesic`](@ref). Its combined return type is
[`KerrGeodesicFamily`](@ref), which can carry stable and plunge branches where
the corresponding branch is available. The older names `KerrGeodesicS` and
`KerrGeodesicSet` remain compatibility aliases.

## Installation

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

## Parameter Conventions

Stable bound orbit calls use APEX-like parameters:

```julia
geo = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
```

Bound plunge calls use constants of motion:

```julia
plunge_family = kerr_geodesic(0.9, (0.94, 0.1, 12.0); radial_start=:turning_point)
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)
```

Do not mix the APEX-like parameter tuple `(a,p,e,x)` with the constants tuple
`(a,E,Lz,Q)`.

## Branch Metadata

Plunge trajectories expose the radial-root class through `OrbitClass` on
[`KerrGeoPlunge`](@ref), and through `RootClass` on [`KerrGeodesicFamily`](@ref).
The implementation labels include `Complex`, `Real1`, and `Real2`. Downstream
code should keep these branch labels in trajectory, source, and waveform
metadata.

Stable-orbit classification metadata is available through
[`kerr_geo_orbit_type_metadata`](@ref). Near the separatrix, inputs within the current
roundoff guard are evaluated at the separatrix radius and labeled `Separatrix`.
Legacy labels do not use `MarginallyStable` or `Unstable`; the structured
metadata still carries machine-facing stability fields.

## Time Coordinates

Bound plunge trajectories expose callable fields:

```julia
t = plunge.Trajectory.t
r = plunge.Trajectory.r
theta = plunge.Trajectory.theta
phi = plunge.Trajectory.phi
rstar = plunge.Trajectory.rstar
u = plunge.Trajectory.u
v = plunge.Trajectory.v
```

The null coordinates follow $u=t-r_*$ and $v=t+r_*$. The optional
`time_origin=:future_horizon_v_zero` policy shifts the coordinate-time origin so
the future-horizon advanced-time anchor is zero for supported bound plunge
plunge trajectories. The retarded time $u$ diverges linearly for an ingoing
future-horizon trajectory; finite plotting coordinates must use an explicitly
recorded cutoff or shifted-display convention.

Bound plunge outputs also include `plunge.Status.duration`, with the
finite Mino-time duration to the event horizon and explicit status fields
recording that Boyer-Lindquist coordinate time and retarded time diverge at the
future horizon.

The near-horizon `u(rstar)` and `v(rstar)` series callables are exposed on the
trajectory object when supported by the branch. The Mino-time trajectory path
remains the owner of `lambda -> r -> rstar`; the near-horizon series is only a
regular advanced-time evaluator and retarded-time proxy at fixed `rstar`.

## Internal Layout

The top-level module keeps the public constructors and structured output types.
Internal plunge helpers are split into `KerrGeoPlunge/NearHorizonTime.jl` for
the tortoise coordinate, near-horizon null-time series, and horizon-anchor
estimation, and `KerrGeoPlunge/InitialConditions.jl` for radial and polar
initial-condition conversion.

## Local Development Note

Project-internal validation scripts in Generic Plunge CodeX use include-style
local loading rather than creating new Julia environments. Package release and
documentation workflows may use standard package-manager steps in GitHub
Actions.
