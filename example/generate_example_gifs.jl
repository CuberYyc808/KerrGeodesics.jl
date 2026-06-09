script_dir = @__DIR__
package_root = normpath(joinpath(script_dir, ".."))
include(joinpath(package_root, "src", "KerrGeodesics.jl"))

using .KerrGeodesics
using GLMakie
using GeometryBasics

function cartesian_from_spherical(r_values, theta_values, phi_values)
    x = r_values .* sin.(theta_values) .* cos.(phi_values)
    y = r_values .* sin.(theta_values) .* sin.(phi_values)
    z = r_values .* cos.(theta_values)
    return x, y, z
end

function add_black_hole!(ax, a; alpha=0.6)
    rplus = 1 + sqrt(1 - a^2)
    sphere_mesh = Sphere(Point3f(0, 0, 0), Float32(rplus))
    mesh!(ax, sphere_mesh; color=:black, transparency=true, alpha=alpha)
    return rplus
end

function set_equal_limits!(ax, x, y, z; pad=0.1)
    maxabs = maximum(abs, vcat(x, y, z))
    lim = maxabs * (1 + pad)
    xlims!(ax, -lim, lim)
    ylims!(ax, -lim, lim)
    zlims!(ax, -lim, lim)
    return nothing
end

function render_stable_gif()
    orbit = kerr_geo_stable(0.9, 10.0, 0.5, 0.8; initPhases=(0.0, 0.0, 0.0, 0.0))
    lambda_values = range(0.0, 20.0, length=900)
    r_values = [orbit.Trajectory.r(lambda) for lambda in lambda_values]
    theta_values = [orbit.Trajectory.θ(lambda) for lambda in lambda_values]
    phi_values = [orbit.Trajectory.ϕ(lambda) for lambda in lambda_values]
    x, y, z = cartesian_from_spherical(r_values, theta_values, phi_values)

    fig = Figure(size=(760, 620))
    ax = Axis3(fig[1, 1], title="Stable bound Kerr geodesic", aspect=:data)
    add_black_hole!(ax, orbit.OrbitalParameters.a; alpha=0.6)
    lines!(ax, x, y, z; color=:steelblue, linewidth=2)
    particle_position = Observable(Point3f[Point3f(x[1], y[1], z[1])])
    scatter!(ax, particle_position; color=:crimson, markersize=12)
    set_equal_limits!(ax, x, y, z)

    frames = range(1, length(x), length=30 * 45)
    record(fig, joinpath(script_dir, "Trajectory_stable.gif"), eachindex(frames); framerate=45) do i
        idx = clamp(round(Int, frames[i]), 1, length(x))
        particle_position[] = Point3f[Point3f(x[idx], y[idx], z[idx])]
    end
    return orbit
end

function render_plunge_gif()
    orbit = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)
    lambda_end = 0.995 * orbit.Status.duration.mino_time_to_horizon
    lambda_values = range(0.0, lambda_end, length=900)
    r_values = [orbit.Trajectory.r(lambda) for lambda in lambda_values]
    theta_values = [orbit.Trajectory.theta(lambda) for lambda in lambda_values]
    phi_values = [orbit.Trajectory.phi(lambda) for lambda in lambda_values]
    x, y, z = cartesian_from_spherical(r_values, theta_values, phi_values)

    fig = Figure(size=(760, 620))
    ax = Axis3(fig[1, 1], title="Bound plunge Kerr geodesic", aspect=:data)
    add_black_hole!(ax, orbit.OrbitalParameters.a; alpha=0.6)
    lines!(ax, x, y, z; color=:darkred, linewidth=2)
    scatter!(ax, [x[1]], [y[1]], [z[1]]; color=:royalblue, markersize=10)
    particle_position = Observable(Point3f[Point3f(x[1], y[1], z[1])])
    scatter!(ax, particle_position; color=:crimson, markersize=12)
    set_equal_limits!(ax, x, y, z)

    frames = range(1, length(x), length=10 * 45)
    record(fig, joinpath(script_dir, "Trajectory_plunge.gif"), eachindex(frames); framerate=45) do i
        idx = clamp(round(Int, frames[i]), 1, length(x))
        particle_position[] = Point3f[Point3f(x[idx], y[idx], z[idx])]
    end
    return orbit
end

function render_example_gifs()
    stable_orbit = render_stable_gif()
    plunge_orbit = render_plunge_gif()
    println("Wrote ", joinpath(script_dir, "Trajectory_stable.gif"))
    println("Wrote ", joinpath(script_dir, "Trajectory_plunge.gif"))
    println("Stable starts at r0 = ", stable_orbit.Trajectory.r(0.0))
    println("Plunge starts at r0 = ", plunge_orbit.Trajectory.r(0.0))
    return stable_orbit, plunge_orbit
end

if abspath(PROGRAM_FILE) == @__FILE__
    render_example_gifs()
end
