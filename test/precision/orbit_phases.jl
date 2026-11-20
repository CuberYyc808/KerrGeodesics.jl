# Mino-time averages against quadrature of their definitions.

const AVERAGE_ORBITS = [(0.9, 6.0, 0.7, 0.5), (0.998, 5.281, 0.0, 4.8e-4), (0.7, 12.0, 0.5, -0.6),
    (0.5, 9.0, 0.0, 0.0), (0.0, 10.0, 0.2, 0.5), (-0.6, 10.0, 0.4, 0.3)]

# mean over one period 2π of a smooth periodic f(q), by the trapezoid rule doubled to convergence
function periodic_mean(f, ::Type{T}, tol) where {T}
    N = 16
    total = sum(f(2 * T(π) * k / N) for k in 0:N-1)
    mean = total / N
    while true
        total += sum(f(2 * T(π) * (2k + 1) / (2N)) for k in 0:N-1)
        N *= 2
        next = total / N
        abs(next - mean) <= tol * max(abs(next), 1) && return next
        N >= 2^19 && error("no convergence")
        mean = next
    end
end

# The averages from the trajectory of kerr_geo_orbit at `bits` of precision: radial ones over
# q_r, polar ones over q_θ.
function reference_averages(c, bits, tol)
    setprecision(BigFloat, bits) do
        T = BigFloat
        o = kerr_geo_orbit(c...; precision=bits)
        a, E, Lz = T(c[1]), o["Energy"], o["AngularMomentum"]
        F = o["Frequencies"]
        r_of(q) = o["Trajectory"][2](q / F["ϒr"])
        θ_of(q) = o["Trajectory"][3](q / F["ϒθ"])
        Δ(r) = r^2 - 2r + a^2
        radial(f) = iszero(F["ϒr"]) ? f(r_of(zero(T))) : periodic_mean(q -> f(r_of(q)), T, tol)
        polar(f) = periodic_mean(q -> f(θ_of(q)), T, tol)
        (z2_theta = polar(θ -> cos(θ)^2),
            lz_over_sin2_theta = iszero(Lz) ? T(NaN) : polar(θ -> Lz / sin(θ)^2),
            phi_r = radial(r -> a * (E * (r^2 + a^2) - a * Lz) / Δ(r)),
            t_r = radial(r -> (r^2 + a^2) * (E * (r^2 + a^2) - a * Lz) / Δ(r)),
            t_theta = polar(θ -> -a^2 * E * sin(θ)^2),
            r2 = radial(r -> r^2))
    end
end

relerr(v, ref) = Float64(abs(v - ref) / max(abs(ref), 1))

@testset "Mino-time averages" begin
    for c in AVERAGE_ORBITS
        m = kerr_geo_mino_averages(c...)
        @test m isa NamedTuple && all(v -> v isa Float64, values(m))
        ref = reference_averages(c, 128, 1e-22)
        for k in keys(m)
            isnan(ref[k]) && continue
            @test relerr(m[k], ref[k]) <= 1e-12
        end
        a = c[1]
        f = kerr_geo_frequencies(c...)
        k = kerr_geo_constants_of_motion(c...)
        @test abs(f["ϒϕ"] - (m.lz_over_sin2_theta + m.phi_r - a * k["E"])) <= 1e-13 * abs(f["ϒϕ"])
        @test abs(f["ϒt"] - (m.t_r + m.t_theta + a * k["Lz"])) <= 1e-13 * f["ϒt"]
        @test abs(KerrGeodesics.kerr_geo_proper_frequency_factor(c...) - (m.r2 + a^2 * m.z2_theta)) <=
            1e-13 * m.r2
    end
    # L_z = 0: the polar average is the limit x → 0± (the orbit passes over the axis)
    a, p, e = 0.5, 9.0, 0.0
    for x in (0.0, -0.0)
        m = kerr_geo_mino_averages(a, p, e, x)
        f = kerr_geo_frequencies(a, p, e, x)
        @test m.lz_over_sin2_theta == (signbit(x) ? -1 : 1) * f["ϒθ"]
    end
    # BigFloat
    c = (0.9, 6.0, 0.7, 0.5)
    m = kerr_geo_mino_averages(c...; precision=256)
    @test all(v -> v isa BigFloat && precision(v) == 256, values(m))
    ref = reference_averages(c, 320, 1e-72)
    for k in keys(m)
        @test relerr(m[k], ref[k]) <= 1e-60
    end
    @test_throws DomainError kerr_geo_mino_averages(0.9, 10.0, 1.0, 0.5)
end
