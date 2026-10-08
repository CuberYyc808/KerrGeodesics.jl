# Concurrent evaluation: different orbits evaluated from several threads at once give the same
# bits as serial evaluation (no shared mutable state; lazily built tables are published
# atomically). The interfaces and one member of every catalogue case are covered.

const LAMBDAS = collect(range(-2.3, 3.1; length=40))

# The evaluations of one orbit: a vector of closures returning numbers, built afresh each call
# so that lazily built tables start empty.
function apex_jobs(c)
    o = kerr_geo_orbit(c...; initPhases=(0.4, 1.1, 2.3, -0.7))
    fs = [o["Trajectory"]; o["FourVelocity"]]
    o["CrossFunction"] === nothing || append!(fs, o["CrossFunction"])
    return fs
end
plunge_jobs(c) = [generic_plunge_orbit(c...; initPhases=(0.2, 0.3, 0.5, -0.1));
                  generic_plunge_velocity(c...; initPhase=(0.3, 0.5))]

const CATALOGUE = [split(l, '\t') for l in eachline(joinpath(@__DIR__, "catalogue_constants.tsv"))]
function member_of(row)
    id, slot, index, a, E, Lz, Q, kwargs = row
    family = kerr_geodesic(parse(Float64, a), (parse(Float64, E), parse(Float64, Lz),
        parse(Float64, Q)); eval(Meta.parse(kwargs))...)
    m = getfield(family, kerr_geo_class(Symbol(slot)).slot)
    return m isa Tuple ? m[parse(Int, index)] : m
end
function member_jobs(row)
    m = member_of(row)
    lo, hi = m.Domain.mino
    lo, hi = max(lo, -5.0), min(hi, 5.0)
    names = filter(n -> n in (:t, :r, :z, :phi, :tau, :v, :psi), propertynames(m.Trajectory))
    return [λ -> getproperty(m.Trajectory, n)(lo + (hi - lo) * (λ + 2.3) / 5.4) for n in names]
end

# values of all closures at all λ; an exception is recorded by its type
function evaluate(jobs)
    out = Any[]
    for f in jobs, λ in LAMBDAS
        push!(out, try f(λ) catch err typeof(err) end)
    end
    return out
end
# bit for bit: === for Float64, the same value and precision for BigFloat (objects)
_same(p, q) = p === q || (p isa BigFloat && q isa BigFloat && precision(p) == precision(q) &&
    isequal(p, q)) || (p isa Number && q isa Number && isnan(p) && isnan(q))
same(x, y) = length(x) == length(y) && all(((p, q),) -> _same(p, q), zip(x, y))

function concurrent_matches_serial(make, cases)
    serial = [evaluate(make(c)) for c in cases]
    for _ in 1:2
        tasks = [Threads.@spawn evaluate(make(c)) for c in cases]
        all(same.(fetch.(tasks), serial)) || return false
    end
    return true
end

# members built at 256 bits (one per class) and the interfaces at 256 bits
function big_member_jobs(row)
    m = setprecision(() -> member_of_big(row), BigFloat, 256)
    lo, hi = Float64.(m.Domain.mino)
    lo, hi = max(lo, -5.0), min(hi, 5.0)
    names = filter(n -> n in (:t, :r, :z, :phi, :tau), propertynames(m.Trajectory))
    return [λ -> getproperty(m.Trajectory, n)(big(lo + (hi - lo) * (λ + 2.3) / 5.4)) for n in names]
end
function member_of_big(row)
    id, slot, index, a, E, L, Q, kwargs = row
    family = kerr_geodesic(big(parse(Float64, a)), (big(parse(Float64, E)), big(parse(Float64, L)),
        big(parse(Float64, Q))); eval(Meta.parse(kwargs))...)
    m = getfield(family, kerr_geo_class(Symbol(slot)).slot)
    return m isa Tuple ? m[parse(Int, index)] : m
end
big_apex_jobs(c) = (o = kerr_geo_orbit(c...; initPhases=(0.4, 1.1, 2.3, -0.7), precision=256);
    [λ -> f(big(λ)) for f in [o["Trajectory"]; o["FourVelocity"]]])

@testset "Concurrent evaluation ($(Threads.nthreads()) threads)" begin
    apex = [(0.9, 10.0, 0.5, 0.8), (0.5, 12.0, 0.3, -0.5), (0.99, 6.0, 0.1, 0.95), (0.0, 12.0, 0.4, 0.6)]
    @test concurrent_matches_serial(apex_jobs, apex)
    plunge = [(0.9, 0.94, 2.0, 4.0), (0.5, 0.97, 1.0, 1.0), (0.9, 0.9, 2.0, 3.0), (0.0, 0.94, 3.0, 1.0)]
    @test concurrent_matches_serial(plunge_jobs, plunge)
    for group in Iterators.partition(CATALOGUE, 4)
        @test concurrent_matches_serial(member_jobs, collect(group))
    end
    big_rows = filter(row -> row[1] in ("A1", "B4", "C3", "D2"), CATALOGUE)
    @test concurrent_matches_serial(big_member_jobs, big_rows)
    @test concurrent_matches_serial(big_apex_jobs, apex[1:2])
end
