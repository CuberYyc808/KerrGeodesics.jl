# KerrGeodesics.jl

A Julia package for computing and visualizing geodesic motion around a Kerr black hole.

## Installation

You can install it directly with HTTPS:

```julia
using Pkg
Pkg.add(url="https://github.com/CuberYyc808/KerrGeodesics.jl.git")
```

---

## Usage

### Basic Example

The main function is `Kerr_Geodesics(a, p, e, x)`, which computes the trajectory of a test particle around a Kerr black hole.

- `a`: Spin parameter of the black hole.
- `p`: Semi-latus rectum of the orbit.
- `e`: Orbital eccentricity.
- `x`: Inclination parameter.

It returns a dictionary containing the trajectory information (`t, r, θ, φ`) as well as the Cartesian coordinates (`x, y, z`).

```julia
using KerrGeodesics

# Example: generic Kerr geodesics with a=0.9, p=10.0, e=0.5, x=0.8. 
assoc = KerrGeoOrbit(0.9, 10.0, 0.5, 0.8; initPhases=(0.0, 0.0, 0.0, 0.0))
```

It should give you a big dictionary

---

## Visualization Example

example/Trajectory_generic.gif

You can find a example of how to visualize your results in [example](example/Test_KerrGeodesics.ipynb)
