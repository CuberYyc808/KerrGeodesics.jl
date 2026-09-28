# Precompile workload. Runs only while the package image is generated, so the machine
# code of the common constructors and trajectory closures is cached: without it the first
# kerr_geodesic call in a session spends ~15 s compiling.
function _precompile_workload()
    constants = (
        (0.9, (0.95, 2.5, 6.0)),                              # plunge (B4)
        (0.7, (1.2, 3.0, 10.0)),                              # hyperbolic plunge + scatter (B6, D2)
        (0.3, (1.0, 4.0, 5.0)),                               # parabolic
        (0.7, (1.2, 0.7176768865182188, 24.938291301882202)), # repeated root (K11/K10)
        (0.0, (0.9622504486493763, 3.6742346141747673, 0.0)), # unstable circular family
        (0.9, (-0.8, -4.0, 1.0)),                             # N4 trapped
        (1.0, (1.2, 2.4, 1.0)),                               # exact extremal
    )
    for (a, c) in constants
        try
            family = kerr_geodesic(a, c)
            for m in kerr_geo_members(family)
                lo, hi = _diagnostic_window(m)
                λ = (lo + hi) / 2
                # property access as users write it (m.Trajectory.t): getproperty is compiled
                # per member NamedTuple type, and these types are large
                for (fields, keys) in ((m.Trajectory, (:r, :t, :phi, :z, :theta, :tau, :v, :psi)),
                        (m.Velocity, (:ut, :ur, :uz, :uphi)))
                    for k in keys
                        hasproperty(fields, k) || continue
                        try getproperty(fields, k)(λ) catch end
                    end
                end
            end
        catch
        end
    end
    # the finite-window (E >= 1) and plunge APIs (their own engine closures)
    for (a, E, L, Q) in ((0.5, 1.1, 0.2, 0.0), (0.9, 1.1, 0.5, 3.0), (0.5, 1.0, 1.0, 0.0),
            (0.7, 1.0, 1.2, 2.0))
        try
            o = kerr_geo_capture(a, E, L, Q; input=:constants)
            rp = o.ReferenceZero.rplus
            λ = 0.5 * o.ReferenceZero.lambda_infinity
            foreach(f -> f(λ), (o.Trajectory.tau, o.Trajectory.v, o.Trajectory.psi))
            foreach(f -> f(rp + 0.2, rp + 0.5), (o.Trajectory.total_time_increment,
                o.Trajectory.total_phi_increment, o.Trajectory.total_v_increment,
                o.Trajectory.total_psi_increment, o.Trajectory.radial_proper_increment))
        catch
        end
    end
    for (a, E, L, Q) in ((0.9, 1.1, 5.0, 3.0), (0.5, 1.1, 5.0, 0.0), (0.5, 1.0, 4.5, 0.0),
            (0.7, 1.0, 4.5, 2.0))
        try
            o = kerr_geo_scatter(a, E, L, Q; input=:constants)
            o.Trajectory.t(0.1); o.Trajectory.phi(0.1)
        catch
        end
    end
    for start in (:turning_point, :inner_turning)
        try
            o = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=start)
            o.Trajectory.t(0.1); o.Trajectory.phi(0.1)
        catch
        end
    end
    try
        kerr_geodesic(0.9, 8.0, 0.3, 0.6)
        kerr_geo_frequencies(0.9, 8.0, 0.3, 0.6; Time="Mino")
        o = kerr_geo_orbit(0.9, 8.0, 0.3, 0.6)
        foreach(f -> f(1.0), o["Trajectory"])
    catch
    end
    return nothing
end

ccall(:jl_generating_output, Cint, ()) == 1 && _precompile_workload()
