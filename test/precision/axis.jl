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

# The jump of φ across one pass over the axis at λc, and the jump of the Cartesian velocity there
# relative to its size (both from one-sided differences at distance δ)
function axis_pass(m, λc; δ=1e-6)
    X(λ) = (r = m.Trajectory.r(λ); θ = m.Trajectory.theta(λ); φ = m.Trajectory.phi(λ);
        (r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)))
    jump = m.Trajectory.phi(λc + δ) - m.Trajectory.phi(λc - δ)
    before = (X(λc - δ) .- X(λc - 2δ)) ./ δ
    after = (X(λc + 2δ) .- X(λc + δ)) ./ δ
    return jump, maximum(abs.(after .- before)) / maximum(abs.(before))
end

# the Mino times in (lo, hi) where sin θ vanishes: the minima of sin²θ on a grid, refined
function axis_passes(m, lo, hi; n=2000)
    s2(λ) = sin(m.Trajectory.theta(λ))^2
    λs = range(lo, hi; length=n)
    v = s2.(λs)
    passes = Float64[]
    for i in 2:n-1
        v[i] <= v[i-1] && v[i] < v[i+1] || continue
        a, b = λs[i-1], λs[i+1]
        for _ in 1:200                      # golden-section search for the minimum
            c = b - (b - a) / Base.MathConstants.golden
            d = a + (b - a) / Base.MathConstants.golden
            s2(c) < s2(d) ? (b = d) : (a = c)
        end
        push!(passes, (a + b) / 2)
    end
    return passes
end

@testset "Lz = 0: φ is continuous over the axis in every member class" begin
    # each case: member, Mino window (inside the domain of the member); φ jumps by +π at every
    # pass and the Cartesian path is smooth
    window(m) = (d = m.Domain.mino; lo = max(d[1], -1.0); hi = min(d[2], 3.0);
        (lo + 0.01 * (hi - lo), hi - 0.01 * (hi - lo)))
    cases = (
        ("Stable", kerr_geodesic(0.9, 8.0, 0.1, 0.0).Stable),
        ("Stable a < 0", kerr_geodesic(-0.9, 8.0, 0.1, 0.0).Stable),
        ("Plunge", kerr_geo_plunge_component(0.5, 0.97, 0.0, 12.0)),
        ("Plunge", kerr_geo_plunge_component(0.9, 0.95, 0.0, 3.0)),
        ("Scatter", kerr_geo_scatter_component(0.9, 1.1, 0.0, 40.0)),
        # Q = 0: z = ±sech, one pass at u = 0, here inside the capture window
        ("Capture Q = 0", kerr_geo_capture_component(0.9, 1.1, 0.0, 0.0;
            polar_sector=:axis_crossing, polar_phase=0.3)),
    )
    for (label, m) in cases
        lo, hi = window(m)
        passes = axis_passes(m, lo, hi)
        @test !isempty(passes)
        for λc in passes
            jump, kink = axis_pass(m, λc)
            @test abs(jump - π) <= 1e-4
            @test kink <= 1e-3                  # a missing π reverses the horizontal velocity: O(1)
        end
        # away from the passes φ changes by much less than π between neighbouring samples
        λs = range(lo, hi; length=4000)
        φ = [m.Trajectory.phi(λ) for λ in λs]
        steps = abs.(diff(φ))
        @test count(>(1.0), steps) == length(passes)
    end
end

@testset "Lz = 0: kerr_geo_orbit at x = ±0.0" begin
    # x = −0.0 is the limit x → 0⁻: φ jumps by −π over the axis and ϒϕ includes −ϒθ
    a, p, e = 0.9, 8.0, 0.1
    fp, fm = kerr_geo_frequencies(a, p, e, 0.0), kerr_geo_frequencies(a, p, e, -0.0)
    op, om = kerr_geo_orbit(a, p, e, 0.0), kerr_geo_orbit(a, p, e, -0.0)
    φp, φm = op["Trajectory"][4], om["Trajectory"][4]
    λc = π / fp["ϒθ"]
    @test φp(λc + 1e-7) - φp(λc - 1e-7) ≈ π atol = 1e-5
    @test φm(λc + 1e-7) - φm(λc - 1e-7) ≈ -π atol = 1e-5
    # the polar part of ϒϕ is ±ϒθ
    @test (fp["ϒϕ"] - fm["ϒϕ"]) / 2 ≈ fp["ϒθ"] rtol = 1e-13
    # the orbit at x = −0.0 agrees with the Stable member near the axis on the negative side
    s = kerr_geodesic(a, p, e, -1e-12).Stable
    @test all(λ -> abs(φm(λ) - s.Trajectory.phi(λ)) <= 1e-9, (0.3, 1.0, 2.3))
    # x = +0.0 is unchanged: the same as the Stable member at x = 0
    s0 = kerr_geodesic(a, p, e, 0.0).Stable
    @test all(λ -> abs(φp(λ) - s0.Trajectory.phi(λ)) <= 1e-12, (0.3, 1.0, 2.3))
end
