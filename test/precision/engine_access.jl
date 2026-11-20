# The public coordinate-engine interface: `engine_coordinates`, `horizon_relative_state`,
# `horizon_relative_rstar`. They forward to the internal methods, so the values are those the
# members use, bit for bit.

@testset "Coordinate engine and horizon-relative state" begin
    for (a, E, L, Q) in ((0.9, 0.95, 2.0, 3.0), (0.5, 0.97, 0.0, 12.0), (0.9, 0.95, -1.0, 0.0),
            (0.999, 0.9, 1.5, 2.0))
        m = kerr_geo_plunge_component(a, E, L, Q)
        c = engine_coordinates(m)
        @test c isa EngineCoordinates
        @test (c.a, c.E, c.L, c.Q) == (a, E, L, Q)
        λH = m.Domain.horizon_lambda
        for λ in (0.0, 0.25λH, 0.75λH, (1 - 1e-6)λH, prevfloat(λH))
            s = horizon_relative_state(c, λ)
            @test s !== nothing
            @test s === horizon_relative_state(c.r_of, λ) === KG._radial_state(c.r_of, λ)
            h = s.chart
            r = h.horizon[1] + (h.horizon[2] + s.gap)
            @test r === m.Trajectory.r(λ)
            @test s.velocity === m.Velocity.ur(λ)
            @test horizon_relative_rstar(s) === KG._relative_rstar(s) === m.Trajectory.rstar(λ)
            @test s.gap >= 0
            # the gap keeps its digits where r − r₊ from r does not
            rp = 1 + sqrt(1 - a^2)
            @test abs(s.gap - (r - rp)) <= 4eps(rp)
            λ <= 0.75λH && @test horizon_relative_rstar(s) ≈ kerr_rstar(a, r) rtol = 1e-12 atol = 1e-12
        end
        @test horizon_relative_state(c, λH).gap == 0
    end
    # members whose radius does not reach the future horizon: no horizon-relative state
    @test horizon_relative_state(r -> 2.0, 0.5) === nothing
    @test engine_coordinates(1.0) === nothing
    s = kerr_geo_scatter_component(0.9, 1.1, 0.0, 40.0)
    c = engine_coordinates(s)
    @test c === nothing || horizon_relative_state(c, 0.1) === nothing
    # the public names are exported and documented
    for name in (:EngineCoordinates, :engine_coordinates, :horizon_relative_state,
            :horizon_relative_rstar)
        @test name in names(KerrGeodesics)
        @test Base.Docs.hasdoc(KerrGeodesics, name)
    end
end
