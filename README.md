# KerrGeodesics.jl

Julia interfaces for Kerr geodesic trajectories in units with `G = c = M = 1`.

Current support:

- stable bound orbits from APEX-like parameters `(a,p,e,x)`;
- bound plunge orbits from constants `(a,E,Lz,Q)`, with `E < 1`.

Scattering orbits are planned for future development.

## Installation

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

## Basic Usage

```julia
using KerrGeodesics

stable = kerr_geo_stable(0.9, 10.0, 0.5, 0.8)
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)
family = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
```

Typical printed outputs are:

```julia
KerrGeoStable(
    OrbitalParameters = (a = 0.9, p = 10.0, e = 0.5, x = 0.8),
    ConstantsOfMotion = (E = ..., Lz = ..., Q = ...),
    OrbitalType = ["Bound", "Eccentric", "Stable", "Inclined"],
    Frequencies = (Υt = ..., Υr = ..., Υθ = ..., Υϕ = ...),
    Parametrization = "Mino",
    Trajectory = (t = t(λ), r = r(λ), θ = θ(λ), ϕ = ϕ(λ)),
    InitialPhases = (qt0 = 0.0, qr0 = 0.0, qθ0 = 0.0, qϕ0 = 0.0),
)
```

```julia
KerrGeoPlunge(
    ConstantsOfMotion = (E = 0.94, Lz = 0.1, Q = 12.0),
    OrbitClass = "Complex",
    Parametrization = "Mino",
    InitialPosition = (t0 = 0.0, r0 = ..., theta0 = ..., phi0 = 0.0),
    Status = (supported = true, reason = "ok", ...),
    Trajectory = (t = t(lambda), r = r(lambda), theta = theta(lambda), phi = phi(lambda), rstar = rstar(lambda), u = u(lambda), v = v(lambda), u_rstar_series = u(rstar), v_rstar_series = v(rstar)),
    Velocity = (ut = ut(lambda), ur = ur(lambda), uz = dz/dlambda, utheta = dtheta/dlambda, uphi = uphi(lambda)),
)
```

```julia
KerrGeodesicFamily(
    InputType = :apex,
    Parameters = (a = 0.9, p = 10.0, e = 0.5, x = 0.8),
    ConstantsOfMotion = (E = ..., Lz = ..., Q = ...),
    RootClass = "Complex",
    HasStable = true,
    HasPlunge = true,
    Status = (stable_supported = true, plunge_supported = true),
)
```

For a stable eccentric orbit, `initPhases=(0,0,0,0)` starts the radial motion at
periapsis, `r(0)=p/(1+e)`. For a bound plunge, `radial_start=:turning_point`
starts at the exterior turning point.

## Examples

The example notebook builds both stable and bound-plunge trajectory animations:

- [`example/Test_KerrGeodesics.ipynb`](example/Test_KerrGeodesics.ipynb)
- [`example/generate_example_gifs.jl`](example/generate_example_gifs.jl)

Generated example images:

![Stable bound Kerr geodesic](example/Trajectory_stable.gif)

![Bound plunge Kerr geodesic](example/Trajectory_plunge.gif)

## Documentation

Documenter.jl sources live under [`docs/`](docs/).

## License

MIT.
