# KerrGeodesics.jl

![license](https://img.shields.io/badge/license-MIT-blue.svg)
[![GitHub release](https://img.shields.io/github/v/release/CuberYyc808/KerrGeodesics.jl.svg)](https://github.com/CuberYyc808/KerrGeodesics.jl/releases)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://CuberYyc808.github.io/KerrGeodesics.jl)

Timelike geodesics outside a Kerr black hole, as functions of Mino time `λ`
(`G = c = M = 1`, Boyer–Lindquist coordinates, `dτ/dλ = Σ = r² + a² cos²θ`).
Give the spin `a` and the constants of motion `(E, Lz, Q)`: the package classifies the radial
motion outside the outer horizon and returns a member for each radial range those constants
allow.

<p align="center">
  <img src="example/animations/showcase_all.gif" width="100%" alt="the 56 orbit cases">
</p>

## Installation

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

## Usage

An orbit is fixed by the spin `a` and the constants of motion `(E, Lz, Q)`, the main input
of `kerr_geodesic`. An orbit with a periapsis can also be given by its APEX parameters
`(a, p, e, x)`. Both inputs are four numbers, so the constants go in a tuple:
`kerr_geodesic(a, (E, Lz, Q))` reads constants, `kerr_geodesic(a, p, e, x)` reads APEX
parameters.

```julia
using KerrGeodesics

family = kerr_geodesic(0.9, (0.9641, 2.8359, 4.5444))   # a, (E, Lz, Q)
family = kerr_geodesic(0.9, 10.0, 0.5, 0.8)             # a, p, e, x: the same orbit

family.Status.case_ids       # (:A1, :B1): a stable orbit and a plunge share these constants
m = family.Stable            # also family.Critical, .Plunge, .Capture, .Scatter, .Trapped
kerr_geo_members(family)     # all of them
```

**Trajectory.** Each coordinate is a function of `λ` on `m.Domain.mino`:

```julia
λ = 1.0
m.Trajectory.t(λ); m.Trajectory.r(λ); m.Trajectory.theta(λ); m.Trajectory.phi(λ); m.Trajectory.tau(λ)
```

Orbits that reach a horizon also give the horizon-regular coordinates `v` and `psi`, finite on
the horizon:

```julia
b = family.Plunge            # from the turning point (λ = 0) into the horizon
b.Trajectory.v(b.Domain.mino[2]), b.Trajectory.psi(b.Domain.mino[2])
```

**Four-velocity.** `m.Velocity` holds the Mino-time rates `ut = dt/dλ`, `ur`, `utheta`,
`uphi` and `dtau_dlambda = Σ`; divide by `dtau_dlambda` for `dx/dτ`:

```julia
m.Velocity.ut(λ) / m.Velocity.dtau_dlambda(λ)   # dt/dτ
```

**Orbital parameters.**

```julia
m.CaseId                 # :A1 (the case; kerr_geo_member_class(m) gives the class, :stable)
m.ConstantsOfMotion      # (a, E, Lz, Q)
m.Roots                  # radial roots and the polar motion
m.Status.apex            # (a, p, e, x), for stable orbits
```

**Frequencies.** Mino-time frequencies of a stable orbit, and the Boyer–Lindquist or proper-time
ones from `(a, p, e, x)`:

```julia
m.Status.frequencies                                          # (ϒt, ϒr, ϒθ, ϒϕ)
kerr_geo_frequencies(0.9, 10.0, 0.5, 0.8; Time="Mino")        # also "BoyerLindquist", "Proper"
```

**Many points at once.**

```julia
s = kerr_geo_sample(m, range(0, 20; length=10_000))   # s.t, s.r, s.theta, s.phi, s.tau, s.ut, …
```

## Examples

[`example/KerrGeodesics_Tutorial.ipynb`](example/KerrGeodesics_Tutorial.ipynb) walks through
the conventions, the classification, one example per class (Stable, Critical, Plunge,
Capture, Scatter, and Trapped with `E < 0`), polar options, the horizon-regular coordinates,
the four-velocity and self-checks (`kerr_geo_diagnose`), the APEX interface and the accuracy
limits. The 56 cases above, with their constants, are in
[`example/KerrGeodesics_56_Orbit_Catalog.ipynb`](example/KerrGeodesics_56_Orbit_Catalog.ipynb).
