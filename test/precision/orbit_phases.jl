# Mino-time averages against quadrature of their definitions; the orbit evaluated by its phases
# against the Dict of kerr_geo_orbit (in Float64 and in BigFloat); finite values at x = 0;
# kerr_geo_orbit_supported against kerr_geo_orbit succeeding; kerr_geo_orbit_shape against
# kerr_geo_orbit_type_metadata.

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

const PHASE_ORBITS = [AVERAGE_ORBITS; (0.9, 3.9, 0.95, 1.0); (0.9, 12.0, 0.5, -1.0); (0.9, 10.0, 0.0, 1.0)]
const PHASES = [2π * k / 32 for k in 0:32]

# the Dict values at the phases q (radial and polar), the velocities from the four-velocity × Σ
function dict_values(o, q)
    F = o["Frequencies"]
    a = o["a"]
    traj, u = o["Trajectory"], o["FourVelocity"]
    cross = o["CrossFunction"] === nothing ? fill(_ -> 0.0, 4) : o["CrossFunction"]
    dcross = o["DerivativesCrossFunction"] === nothing ? fill(_ -> 0.0, 4) : o["DerivativesCrossFunction"]
    Σ(λ) = traj[2](λ)^2 + a^2 * cos(traj[3](λ))^2
    λr, λθ = q / F["ϒr"], q / F["ϒθ"]
    radial = (r = iszero(F["ϒr"]) ? traj[2](0) : traj[2](λr),
        dr_dlambda = iszero(F["ϒr"]) ? zero(q) : u[2](λr) * Σ(λr),
        delta_t = cross[1](q), delta_phi = cross[3](q), ddelta_t_dq = dcross[1](q),
        ddelta_phi_dq = dcross[3](q))
    θ = traj[3](λθ)
    polar = (z = cos(θ), theta = θ, dtheta_dlambda = u[3](λθ) * Σ(λθ), delta_t = cross[2](q),
        delta_phi = cross[4](q), ddelta_t_dq = dcross[2](q), ddelta_phi_dq = dcross[4](q))
    return radial, polar
end

# largest |value − reference| over the phases, relative to the largest |reference| (or 1)
function phase_error(values, refs, k)
    v = [getfield(x, k) for x in values]
    r = [getfield(x, k) for x in refs]
    good = isfinite.(r)
    any(good) || return 0.0
    return Float64(maximum(abs.(v[good] .- r[good])) / max(maximum(abs.(r[good])), 1))
end

@testset "Orbit evaluated by its phases" begin
    for c in PHASE_ORBITS
        o = kerr_geo_orbit_phases(c...)
        @test o isa KerrGeoOrbitPhases{Float64} && isconcretetype(typeof(o))
        D = kerr_geo_orbit(c...)
        @test (o.E, o.Lz, o.Q) == (D["Energy"], D["AngularMomentum"], D["CarterConstant"])
        @test (o.Upsilon_r, o.Upsilon_theta, o.Upsilon_phi, o.Upsilon_t) ==
            Tuple(D["Frequencies"][k] for k in ("ϒr", "ϒθ", "ϒϕ", "ϒt"))
        rad = [kerr_geo_radial_motion(o, q) for q in PHASES]
        pol = [kerr_geo_polar_motion(o, q) for q in PHASES]
        @test eltype(rad) <: NamedTuple && all(x -> all(v -> v isa Float64 && isfinite(v), values(x)), rad)
        @test all(x -> all(v -> v isa Float64 && isfinite(v), values(x)), pol)
        # the Dict in BigFloat: the record has the accuracy of Float64
        B = kerr_geo_orbit(c...; precision=128)
        refs = setprecision(() -> [dict_values(B, big(q)) for q in PHASES], BigFloat, 128)
        for k in keys(refs[1][1])
            @test phase_error(rad, first.(refs), k) <= 1e-12
        end
        # at L_z = 0 the polar azimuth steps on the passes over the axis, q = jπ: Float64(π) lies
        # on the pass in Float64 and before it in BigFloat: those phases are compared in Float64 below
        off_axis = [!(iszero(o.Lz) && (i - 1) % 16 == 0) for i in eachindex(PHASES)]
        for k in keys(refs[1][2])
            sel = k == :delta_phi ? off_axis : trues(length(PHASES))
            @test phase_error(pol[sel], last.(refs)[sel], k) <= 1e-12
        end
        # the Dict in Float64: the same values within its own rounding errors, which grow as
        # the polar characteristic 1 − z₋² and sin θ (from z) lose digits on a nearly polar orbit
        refs64 = [dict_values(D, q) for q in PHASES]
        refsB = setprecision(() -> [dict_values(B, big(q)) for q in PHASES], BigFloat, 128)
        for (i, vals) in ((1, rad), (2, pol))
            r64, rB = getindex.(refs64, i), getindex.(refsB, i)
            for k in keys(r64[1])
                @test phase_error(vals, r64, k) <= max(1e-12, 2 * phase_error(r64, rB, k))
            end
        end
        # dr/dλ = ±√R(r), dz/dλ = ±√Θ(z)
        a, E, Lz, Q = c[1], o.E, o.Lz, o.Q
        R(r) = (E * (r^2 + a^2) - a * Lz)^2 - (r^2 - 2r + a^2) * (r^2 + (Lz - a * E)^2 + Q)
        Θ(z) = Q - z^2 * (a^2 * (1 - E^2) * (1 - z^2) + Lz^2 + Q)
        scaleR = maximum(abs(R(x.r)) for x in rad) + maximum(x.dr_dlambda^2 for x in rad) + eps()
        scaleΘ = Q + eps()
        if !iszero(o.Upsilon_r)
            @test maximum(abs(x.dr_dlambda^2 - R(x.r)) for x in rad) <= 1e-12 * scaleR
        end
        @test maximum(abs(x.dz_dlambda^2 - Θ(x.z)) for x in pol) <= 1e-12 * scaleΘ
        # exact limits: e = 0 constant radius, |x| = 1 equatorial
        if c[3] == 0
            @test all(x -> x.r == c[2] && x.dr_dlambda == 0 && x.delta_t == 0 && x.delta_phi == 0 &&
                x.ddelta_t_dq == 0 && x.ddelta_phi_dq == 0, rad)
        end
        if abs(c[4]) == 1
            @test all(x -> x.z == 0 && x.theta == π / 2 && x.dz_dlambda == 0 && x.dtheta_dlambda == 0 &&
                x.delta_t == 0 && x.delta_phi == 0 && x.ddelta_t_dq == 0 && x.ddelta_phi_dq == 0, pol)
        end
        # vector methods
        vr = kerr_geo_radial_motion(o, PHASES)
        vp = kerr_geo_polar_motion(o, PHASES)
        @test all(k -> getfield(vr, k) == [getfield(x, k) for x in rad], keys(vr))
        @test all(k -> getfield(vp, k) == [getfield(x, k) for x in pol], keys(vp))
        @test kerr_geo_radial_motion!(vr, o, reverse(PHASES)).r == reverse([x.r for x in rad])
    end
    @test_throws ArgumentError kerr_geo_orbit_phases(0.9, 2.0, 0.3, 0.5)
    @test_throws ArgumentError kerr_geo_orbit_phases(0.9, 10.0, 1.5, 0.5)
    o = kerr_geo_orbit_phases(0.9, 6.0, 0.7, 0.5)
    @test_throws ArgumentError kerr_geo_radial_motion!((r = zeros(3),), o, zeros(3))
end

@testset "Orbit by its phases at x = 0" begin
    δ = 1e-8
    for (a, p, e) in ((0.9, 8.0, 0.0), (0.9, 8.0, 0.3))
        o0, op, om = (kerr_geo_orbit_phases(a, p, e, x) for x in (0.0, δ, -δ))
        for (f, keys_) in ((kerr_geo_radial_motion, (:r, :dr_dlambda, :delta_t, :delta_phi,
                :ddelta_t_dq, :ddelta_phi_dq)), (kerr_geo_polar_motion, (:z, :dz_dlambda,
                :dtheta_dlambda, :delta_t, :ddelta_t_dq)))
            v0, vp, vm = ([f(o, q) for q in PHASES] for o in (o0, op, om))
            @test all(x -> all(isfinite, values(x)), v0)
            for k in keys_
                scale = max(maximum(abs(getfield(x, k)) for x in v0), 1)
                mean = [(getfield(x, k) + getfield(y, k)) / 2 for (x, y) in zip(vp, vm)]
                # E, L_z, Q change linearly in x, so each side differs by O(x); their mean by O(x²)
                @test maximum(abs.(mean .- [getfield(x, k) for x in v0])) <= 1e-12 * scale
                @test maximum(abs(getfield(x, k) - getfield(y, k)) for (x, y) in zip(vp, v0)) <= 1e-6 * scale
            end
        end
        # the polar azimuth: the side x → 0± away from the passes over the axis, the mean on them
        for (x0, x1) in ((0.0, δ), (-0.0, -δ))
            o0, o1 = kerr_geo_orbit_phases(a, p, e, x0), kerr_geo_orbit_phases(a, p, e, x1)
            for (i, q) in enumerate(PHASES)
                (i - 1) % 16 == 0 && continue
                @test abs(kerr_geo_polar_motion(o0, q).delta_phi - kerr_geo_polar_motion(o1, q).delta_phi) <= 1e-6
            end
            side = signbit(x0) ? -1 : 1
            @test kerr_geo_polar_motion(o0, π).delta_phi == 0
            @test kerr_geo_polar_motion(o0, π).dtheta_dlambda == 0
            @test kerr_geo_polar_motion(o0, nextfloat(Float64(π))).delta_phi ≈ side * π / 2
            @test kerr_geo_polar_motion(o0, prevfloat(Float64(π))).delta_phi ≈ -side * π / 2
        end
    end
end

scalar_allocations(o) = (kerr_geo_radial_motion(o, 0.3); kerr_geo_polar_motion(o, 0.3);
    (@allocated kerr_geo_radial_motion(o, 0.3)) + (@allocated kerr_geo_polar_motion(o, 0.3)))
vector_allocations(o, out, q) = (kerr_geo_polar_motion!(out, o, q); @allocated kerr_geo_polar_motion!(out, o, q))

@testset "Orbit by its phases: allocations, tasks, BigFloat" begin
    for c in ((0.9, 6.0, 0.7, 0.5), (0.5, 9.0, 0.0, 0.0), (0.9, 12.0, 0.5, -1.0))
        o = kerr_geo_orbit_phases(c...)
        @test scalar_allocations(o) == 0
        q = collect(range(-3.0, 9.0; length=1000))
        @test vector_allocations(o, kerr_geo_polar_motion(o, q), q) == 0
        serial = (kerr_geo_radial_motion(o, q), kerr_geo_polar_motion(o, q))
        tasks = [Threads.@spawn((kerr_geo_radial_motion(o, q), kerr_geo_polar_motion(o, q))) for _ in 1:2]
        @test all(t -> fetch(t) == serial, tasks)
    end
    c = (0.9, 6.0, 0.7, 0.5)
    o = kerr_geo_orbit_phases(c...; precision=256)
    @test o isa KerrGeoOrbitPhases{BigFloat} && o.precision == 256
    B = kerr_geo_orbit(c...; precision=320)
    for q in (0.0, 0.4, 2.0, 4.1)
        rad, pol = kerr_geo_radial_motion(o, q), kerr_geo_polar_motion(o, q)
        @test precision(rad.r) == 256 && precision(pol.z) == 256
        ref = setprecision(() -> dict_values(B, big(q)), BigFloat, 320)
        for k in keys(ref[1])
            @test relerr(getfield(rad, k), getfield(ref[1], k)) <= 1e-60
        end
        for k in keys(ref[2])
            @test relerr(getfield(pol, k), getfield(ref[2], k)) <= 1e-60
        end
    end
end

# a grid with x = ±1, the separatrix neighbourhood, circular orbits inside the ISSO, E ≥ 1 and
# parameters outside the domain
function predicate_grid()
    grid = Tuple{Float64,Float64,Float64,Float64}[]
    for a in (0.0, 0.5, 0.9, -0.7, 0.999, 1.0), x in (-1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1 - 1e-13),
            e in (0.0, 1e-13, 0.3, 0.9, 1.0, 1.2), p in (0.5, 1.5, 3.0, 4.0, 6.0, 10.0, 20.0)
        push!(grid, (a, p, e, x))
    end
    for a in (0.0, 0.9, -0.7), x in (-1.0, 0.0, 0.5, 1.0), e in (0.0, 0.3, 0.9)
        ps = try e == 0 ? kerr_geo_isso(a, x) : kerr_geo_separatrix(a, e, x) catch; continue end
        for s in (-1e-6, -1e-12, -1e-15, 0.0, 1e-15, 1e-12, 1e-6)
            push!(grid, (a, ps * (1 + s), e, x))
        end
        if e == 0
            rph = KerrGeodesics.kerr_geo_photon_sphere_radius(a, x)
            ri = kerr_geo_ibso(a, x)
            for t in (0.1, 0.5, 0.9)
                push!(grid, (a, ri + t * (ps - ri), 0.0, x), (a, rph + t * (ri - rph), 0.0, x))
            end
        end
    end
    append!(grid, [(1.1, 10.0, 0.3, 0.5), (0.5, -1.0, 0.3, 0.5), (0.5, 10.0, -0.1, 0.5),
        (0.5, 10.0, -1e-13, 0.5), (0.5, 10.0, 0.3, 1.1), (0.5, 10.0, 0.0, 1 + 1e-13), (0.5, 10.0, 1.5, 0.5)])
    return grid
end

@testset "kerr_geo_orbit_supported" begin
    builds(c) = try kerr_geo_orbit(c...); true catch err; err isa Union{ArgumentError,DomainError} || rethrow(); false end
    nsupported = 0
    for c in predicate_grid()
        s = kerr_geo_orbit_supported(c...)
        @test s isa Bool
        @test s == builds(c)
        nsupported += s
    end
    @test nsupported > 100
    # within the tolerance of e = 0 inside the ISSO: circular for the classification, but only
    # e = 0 is a constant-radius orbit
    @test kerr_geo_orbit_supported(0.5, 4.0, 0.0, 0.5) && !kerr_geo_orbit_supported(0.5, 4.0, 1e-13, 0.5)
    @test_throws ArgumentError kerr_geo_orbit(0.5, 4.0, 1e-13, 0.5)
end

@testset "kerr_geo_orbit_shape" begin
    nshape = 0
    for c in predicate_grid()
        shape = kerr_geo_orbit_shape(c[3], c[4])
        @test shape in (:circular, :eccentric, :inclined, :generic)
        meta = try kerr_geo_orbit_type_metadata(c...) catch; nothing end
        meta === nothing && continue
        @test (meta.shape == "Circular") == (shape in (:circular, :inclined))
        @test (meta.inclination == "Equatorial") == (shape in (:circular, :eccentric))
        nshape += 1
    end
    @test nshape > 1000
    @test kerr_geo_orbit_shape(1e-13, 1 - 1e-13) == :circular
    @test kerr_geo_orbit_shape(2e-12, -1.0) == :eccentric
    @test kerr_geo_orbit_shape(0.0, 0.5) == :inclined
    @test kerr_geo_orbit_shape(0.3, -0.5) == :generic
    @test kerr_geo_orbit_shape(big(1e-13), 1.0) == :eccentric      # tolerance of 256-bit BigFloat
end
