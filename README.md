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

# Example: Schwarzschild case, circular orbit at r=10M
assoc = KerrGeoOrbit(0.0, 10.0, 0.0, 0.0)
println(keys(assoc))
```

---

## Visualization Example

We use [GLMakie.jl](https://makie.juliaplots.org/stable/) to plot and animate the trajectory.

```julia
using KerrGeodesics, GLMakie, Interpolations

assoc = KerrGeoOrbit(0.9, 6.0, 0.2, 0.5)
x_values, y_values, z_values = assoc["x"], assoc["y"], assoc["z"]
t_values = assoc["t"]
a = assoc["a"]
rp = 1 + sqrt(1 - a^2)

nframes = 300
fps = 30

# Interpolation
itp_x = LinearInterpolation(t_values, x_values)
itp_y = LinearInterpolation(t_values, y_values)
itp_z = LinearInterpolation(t_values, z_values)
t_video = range(t_values[1], t_values[end], length=nframes)

# Figure
fig = Figure(size=(800, 600))
ax = Axis3(fig[1, 1], title="Geodesic trajectory around a Kerr black hole")
sphere_mesh = Sphere(Point3f(0,0,0), rp)
mesh!(ax, sphere_mesh, color=:black)

lines!(ax, x_values, y_values, z_values, color=:gray)

point = Observable(Point3f(itp_x(t_video[1]), itp_y(t_video[1]), itp_z(t_video[1])))
scatter!(ax, point; color=:red, markersize=10)

record(fig, "Trajectory_generic.gif", 1:nframes; framerate=fps) do i
    point[] = Point3f(itp_x(t_video[i]), itp_y(t_video[i]), itp_z(t_video[i]))
end
```

The output will be an animated GIF like this:

![Trajectory](Trajectory_generic.gif)

---

## Testing

A basic smoke test is included in `test/runtests.jl`:

```julia
using Test
using KerrGeodesics

@testset "basic" begin
    vals = KerrGeoOrbit(0.0, 10.0, 0.5, 0.8)
    @test isa(vals, Dict) || isa(vals, NamedTuple)
end
```

Run with:

```julia
Pkg.test("KerrGeodesics")
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
