# Orbit quantities in the precision of the constants: the radial increments of catalogue
# members between two radii (Mino time, t, φ, τ) against mpmath quadrature.

const INCREMENT_CASES = [split(l, '\t') for l in eachline(joinpath(@__DIR__, "increment_cases.tsv"))]
const INCREMENT_REFS = Dict(first(f) => f[2:end] for f in
    (split(l, '\t') for l in eachline(joinpath(@__DIR__, "increment_refs.tsv"))))
_bits(s) = startswith(s, "0x") ? reinterpret(Float64, parse(UInt64, s[3:end]; base=16)) : parse(Float64, s)

# the member of a catalogue row built in type T (binary64 constants, exact in every precision)
function catalogue_member(::Type{T}, slot, index, a, E, L, Q, kwargs) where {T}
    family = kerr_geodesic(T(_bits(a)), (T(_bits(E)), T(_bits(L)), T(_bits(Q)));
        eval(Meta.parse(kwargs))...)
    m = getfield(family, kerr_geo_class(Symbol(slot)).slot)
    return m isa Tuple ? m[parse(Int, index)] : m
end

# largest relative error of the four radial increments over the cases, at `bits`
function increment_error(bits)
    worst = 0.0
    T = bits == 53 ? Float64 : BigFloat
    setprecision(BigFloat, max(bits, 64)) do
        for (id, slot, index, a, E, L, Q, kw, r1s, r2s) in INCREMENT_CASES
            m = catalogue_member(T, slot, index, a, E, L, Q, kw)
            r1, r2 = T(_bits(r1s)), T(_bits(r2s))
            tr = m.Trajectory
            values = (tr.radial_mino_increment(r1, r2), tr.radial_time_increment(r1, r2),
                tr.radial_phi_increment(r1, r2), tr.radial_proper_increment(r1, r2))
            @test all(v -> v isa T, values)
            for (v, ref) in zip(values, INCREMENT_REFS[id])
                exact = setprecision(() -> parse(BigFloat, ref), BigFloat, 600)
                err = setprecision(() -> Float64(abs(big(v) - exact) / max(abs(exact), 1)), BigFloat, 600)
                worst = max(worst, err)
            end
        end
    end
    return worst
end

@testset "Radial increments of catalogue members against mpmath" begin
    @test increment_error(53) <= 1e-12
    @test increment_error(128) <= 1e-28
    @test increment_error(256) <= 1e-60
end

include("exact_constants.jl")

# A member of every class, built at 128 and 256 bits from exact constants of that precision:
# its numbers are BigFloat and the two precisions agree to the 128-bit level (no quantity is
# held at Float64 accuracy).
const LADDER_IDS = ("A1", "A_H1", "A_X1", "K4", "K8", "K11", "B3", "B5", "B_X1", "C1", "C3",
    "C5", "C_X2", "D2", "D_X1", "N1")
function ladder_member(bits, row)
    id, slot, index, a, E, L, Q, kw = row
    setprecision(BigFloat, bits) do
        c = exact_catalogue_constants(BigFloat, Symbol(id), parse.(Float64, (a, E, L, Q))...)
        family = kerr_geodesic(c[1], (c[2], c[3], c[4]); eval(Meta.parse(kw))...)
        m = getfield(family, kerr_geo_class(Symbol(slot)).slot)
        m isa Tuple ? m[parse(Int, index)] : m
    end
end

@testset "Precision ladder of members of every class" begin
    for row in (split(l, '\t') for l in eachline(joinpath(@__DIR__, "catalogue_constants.tsv")))
        row[1] in LADDER_IDS || continue
        m128, m256 = ladder_member(128, row), ladder_member(256, row)
        lo, hi = Float64.(m256.Domain.mino)
        lo, hi = max(lo, -2.0), min(hi, 2.0)
        for s in (0.2, 0.6), name in (:t, :r, :z, :phi, :tau)
            λ = lo + (hi - lo) * s
            v128 = setprecision(() -> getproperty(m128.Trajectory, name)(big(λ)), BigFloat, 128)
            v256 = setprecision(() -> getproperty(m256.Trajectory, name)(big(λ)), BigFloat, 256)
            @test v256 isa BigFloat && precision(v256) == 256
            @test abs(v128 - v256) <= 1e-30 * max(1, abs(v256))
        end
        @test isempty(non_bigfloat((m256.ConstantsOfMotion, m256.Domain, m256.Roots,
            m256.ReferenceZero)))
    end
end
