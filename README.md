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

It returns a dictionary containing all the information of the set of the parameters `(a, p, e, x)`.

```julia
using KerrGeodesics

# Example: generic Kerr geodesics with a=0.9, p=10.0, e=0.5, x=0.8. 
assoc = Kerr_Geodesics(0.9, 10.0, 0.5, 0.8; initPhases=(0.0, 0.0, 0.0, 0.0))
```

The output should be like:

`
Dict{String, Any} with 18 entries:
  "RadialFrequency"    => 2.79272
  "e"                  => 0.5
  "Cosθ_inc"           => 0.8
  "a"                  => 0.9
  "AngularMomentum"    => 2.83592
  "Frequencies"        => Dict("ϒϕ"=>3.73576, "ϒt"=>171.093, "ϒr"=>2.79272, "ϒθ"=>3.55291)
  "InitialPhases"      => (0.0, 0.0, 0.0, 0.0)
  "RadialRoots"        => [20.0, 6.66667, 1.44177, 0.271714]
  "Energy"             => 0.96412
  "ConstantsOfMotion"  => Dict("Q"=>4.54441, "Lz"=>2.83592, "E"=>0.96412)
  "PolarFrequency"     => 3.55291
  "CarterConstant"     => 4.54441
  "Parametrization"    => "Mino"
  "FourVelocity"       => Function[ut_contrav, ur_contrav, uθ_contrav, uφ_contrav]
  "Trajectory"         => Function[t, r, θ, ϕ]
  "Type"               => ["Bound", "Eccentric", "Inclined"]
  "AzimuthalFrequency" => 3.73576
  "p"                  => 10.0
`

---

## Visualization Example

![Particle trajectory around Kerr black hole](example/Trajectory_generic.gif)

You can find a example of how to visualize your results in [example](example/Test_KerrGeodesics.ipynb)
