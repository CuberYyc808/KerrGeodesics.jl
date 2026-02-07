# KerrGeodesics.jl

This Julia package is a reimplementation of the Mathematica code [KerrGeodesics](https://github.com/BlackHolePerturbationToolkit/KerrGeodesics), adapted and optimized for Julia.

## Installation

You can install the package by simply typing 

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

---

## Usage

### Basic Example

The main function is `Kerr_Geodesics(a, p, e, x)`, which computes the trajectory of a test particle around a Kerr black hole.

- `a`: Spin parameter of the black hole.
- `p`: Semi-latus rectum of the orbit.
- `e`: Orbital eccentricity.
- `x`: Inclination parameter.

It returns a dictionary containing all the information of the set of the parameters `(a, p, e, x)`.

```julia
using KerrGeodesics

# Example: generic Kerr geodesics with a=0.9, p=10.0, e=0.5, x=0.8 and initial phases (0.0, 0.0, 0.0, 0.0). 
EMRI_info = kerr_geo_emri(0.9, 10.0, 0.2, 0.8)
```

The output should be like:

```julia
KerrGeoEMRI(
    OrbitalParameters = (a = 0.9, p = 10.0, e = 0.2, x = 0.8),
    ConstantsOfMotion = (E = 0.9546178536960606, Lz = 2.810660381278431, Q = 4.469510431717436),
    OrbitalType = ["Bound", "Eccentric", "Inclined"],
    Frequencies = (ϒt = 124.77891227683266, ϒr = 2.7729739639615034, ϒθ = 3.521699666954832, ϒϕ = 3.6986834364993015),
    Parametrization = "Mino",
    Trajectory = (t = t(λ), r = r(λ), θ = θ(λ), ϕ = ϕ(λ)),
    InitialPhases = (qt0 = 0.0, qr0 = 0.0, qθ0 = 0.0, qϕ0 = 0.0),
)
```

---

## Visualization Example

![Particle trajectory around Kerr black hole](example/Trajectory_generic.gif)

You can find an example of how to visualize your results in [example](example/Test_KerrGeodesics.ipynb)

## License
The package is licensed under the MIT License.
