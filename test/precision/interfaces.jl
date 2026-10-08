# The APEX and plunge interfaces in the precision of their inputs: constants of motion and Mino
# frequencies against mpmath (100 digits), the precision keyword, the precision an object
# remembers, and a 128/256-bit ladder of the trajectories.

const APEX_REFS = [split(l, '\t') for l in eachline(joinpath(@__DIR__, "apex_refs.tsv"))]

# largest relative error of E, Lz, Q, ϒr, ϒθ over the cases, at `bits`
function apex_error(bits)
    worst = 0.0
    for row in APEX_REFS
        c = parse.(Float64, row[1:4])
        refs = setprecision(() -> parse.(BigFloat, row[5:9]), BigFloat, 400)
        k, f = setprecision(BigFloat, max(bits, 64)) do
            T = bits == 53 ? Float64 : BigFloat
            kerr_geo_constants_of_motion(T.(c)...), kerr_geo_frequencies(T.(c)...)
        end
        values = (k["E"], k["Lz"], k["Q"], f["ϒr"], f["ϒθ"])
        @test all(v -> v isa (bits == 53 ? Float64 : BigFloat), values)
        for (v, r) in zip(values, refs)
            worst = max(worst, setprecision(() -> Float64(abs(big(v) - r) / max(abs(r), 1)), BigFloat, 400))
        end
    end
    return worst
end

@testset "APEX constants and frequencies against mpmath" begin
    @test apex_error(53) <= 1e-13
    @test apex_error(128) <= 1e-33
    @test apex_error(256) <= 1e-70
end

@testset "Precision keyword and remembered precision" begin
    o = kerr_geodesic(0.9, (0.94, 2.0, 4.0); precision=256)
    v = setprecision(() -> o.Plunge.Trajectory.t(big(0.03)), BigFloat, 64)
    @test v isa BigFloat && precision(v) == 256
    @test o.Plunge.ConstantsOfMotion.E isa BigFloat
    # a rational input is rounded once, in the target precision
    r = kerr_geodesic(9 // 10, (47 // 50, 2, 4); precision=128)
    @test r.Plunge.ConstantsOfMotion.E == setprecision(() -> BigFloat(47 // 50), BigFloat, 128)
    orbit = kerr_geo_orbit(0.9, 10.0, 0.5, 0.8; precision=256)
    w = setprecision(() -> orbit["Trajectory"][1](big(0.4)), BigFloat, 80)
    @test precision(w) == 256
    plunge = kerr_geo_plunge(0.9, 0.94, 2.0, 4.0; precision=256)
    @test precision(setprecision(() -> plunge.Trajectory.r(big(0.02)), BigFloat, 64)) == 256
    @test kerr_geo_frequencies(0.9, 10.0, 0.5, 0.8; precision=256)["ϒr"] isa BigFloat
    @test kerr_geo_constants_of_motion(0.9, 10.0, 0.5, 0.8; precision=256)["E"] isa BigFloat
    # Float64 keeps Float64 objects
    @test kerr_geodesic(0.9, (0.94, 2.0, 4.0)).Plunge.Trajectory.t(0.03) isa Float64
end

@testset "Interface trajectories at 128 and 256 bits" begin
    for c in ((0.9, 10.0, 0.5, 0.8), (0.5, 7.0, 0.3, 0.5))
        o128, o256 = kerr_geo_orbit(c...; precision=128), kerr_geo_orbit(c...; precision=256)
        for (f128, f256) in zip([o128["Trajectory"]; o128["FourVelocity"]], [o256["Trajectory"]; o256["FourVelocity"]]),
                λ in big.((-1.3, 0.7))
            v128, v256 = f128(λ), f256(λ)
            @test abs(v128 - v256) <= 1e-30 * max(1, abs(v256))
        end
    end
    p128, p256 = kerr_geo_plunge(0.9, 0.94, 2.0, 4.0; precision=128),
        kerr_geo_plunge(0.9, 0.94, 2.0, 4.0; precision=256)
    for name in (:t, :r, :theta, :phi), λ in big.((0.01, 0.05))
        v128, v256 = getproperty(p128.Trajectory, name)(λ), getproperty(p256.Trajectory, name)(λ)
        @test abs(v128 - v256) <= 1e-30 * max(1, abs(v256))
    end
end

@testset "Near-horizon series of a Real2 plunge" begin
    c = (0.9, 0.94, 1.2, 0.01)                  # Real2 with polar motion
    exact(x) = setprecision(BigFloat, 256) do
        o = kerr_geo_plunge(big.(c)...; theta_phase=big(0.3))
        a = big(c[1]); rp = 1 + sqrt(1 - a^2); r = rp + big(x) * rp
        λ = lambda_of_r(big.(c)...)[3](r) - o.InitialPhases.radial
        rs = kerr_rstar(a, r)
        (o.Trajectory.t(λ) + rs, o.Trajectory.v_rstar_series(rs), rs)
    end
    o64 = kerr_geo_plunge(c...; theta_phase=0.3)
    for x in (1e-6, 1e-4, 1e-2)
        v, series, rs = exact(x)
        @test abs(series - v) <= 1e-70 * abs(v)
        @test abs(o64.Trajectory.v_rstar_series(Float64(rs)) - v) <= 1e-13 * abs(v)
    end
end
