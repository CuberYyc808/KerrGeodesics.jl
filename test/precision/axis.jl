# Orbits over the spin axis. At Lz = 0 the azimuth is the limit Lz → 0⁺: φ gains π at every pass
# over the axis, so the frequencies include ϒθ and the Cartesian path continues to the opposite
# meridian. For a = 0 the orbit is planar and φ(λ) is known in closed form at every x.

# a = 0: in-plane angle ψ = π/2 + Lλ (L = p/√(p − 3 − e²)), turning point nearest the north pole
# at λ = 0, φ = atan2(x sin ψ, cos ψ) − π/2 continued through the nodes
function planar_phi(p, e, x, λ)
    L = p / sqrt(p - 3 - e^2)
    π_ = oftype(L, π)
    ψ = π_ / 2 + L * λ
    raw = atan(x * sin(ψ), cos(ψ))
    return raw + 2π_ * round((ψ - raw) / (2π_)) - π_ / 2
end

stable(a, p, e, x; bits=nothing) = bits === nothing ?
    kerr_geodesic(a, p, e, x; initPhases=(0, 0, 0, 0)).Stable :
    kerr_geodesic(a, p, e, x; precision=bits, initPhases=(0, 0, 0, 0)).Stable

@testset "Lz = 0: φ is the limit Lz → 0⁺" begin
    # a = 0 against the planar orbit, at x = 0 and in the near-axis range Lz ≲ 1e-18 at 128 bits
    for x in (0, big(1) // big(10)^25, big(1) // 2)
        m = stable(0 // 1, 7 // 1, 1 // 10, x; bits=128)
        setprecision(BigFloat, 192) do
            for λ in (1 // 3, 1 // 1, 7 // 3)
                v = setprecision(() -> m.Trajectory.phi(BigFloat(λ)), BigFloat, 128)
                ref = planar_phi(big(7), big(1) / 10, big(x), big(λ))
                @test abs(v - ref) / (abs(ref) + 1) <= 1e-36
            end
        end
    end
    m = stable(0.0, 7.0, 0.0, 0.0)
    @test m.Status.frequencies.ϒϕ ≈ m.Status.frequencies.ϒθ rtol = 4eps()
    @test m.Trajectory.phi(1.0) ≈ 3π / 2 rtol = 4eps()
    for a in (0.9, -0.9)
        m0, m1 = stable(a, 8.0, 0.1, 0.0), stable(a, 8.0, 0.1, 1e-12)
        @test m0.Status.frequencies.ϒϕ ≈ m1.Status.frequencies.ϒϕ rtol = 1e-10
        @test all(λ -> abs(m0.Trajectory.phi(λ) - m1.Trajectory.phi(λ)) <= 1e-10, (0.3, 1.0, 2.3))
        f = kerr_geo_frequencies(a, 8.0, 0.1, 0.0)
        @test f["ϒϕ"] ≈ m0.Status.frequencies.ϒϕ rtol = 1e-13
        @test f["ϒϕ"] ≈ kerr_geo_frequencies(a, 8.0, 0.1, 1e-12)["ϒϕ"] rtol = 1e-10
        orbit = kerr_geo_orbit(a, 8.0, 0.1, 0.0)
        @test all(λ -> abs(orbit["Trajectory"][4](λ) - m0.Trajectory.phi(λ)) <= 1e-12, (0.3, 1.0, 2.3))
        # Cartesian velocity is continuous across the first pass over the south pole
        λc = π / m0.Status.frequencies.ϒθ
        X(λ) = (r = m0.Trajectory.r(λ); θ = m0.Trajectory.theta(λ); φ = m0.Trajectory.phi(λ);
            (r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)))
        δ = 1e-6
        before = (X(λc - δ) .- X(λc - 2δ)) ./ δ
        after = (X(λc + 2δ) .- X(λc + δ)) ./ δ
        @test maximum(abs.(after .- before)) <= 1e-4 * maximum(abs.(before))
    end
end
