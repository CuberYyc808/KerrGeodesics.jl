# Working precision: the tolerance and truncation helpers leave Float64 unchanged; the elliptic
# integrals and Jacobi functions reach the precision of their arguments (mpmath references at
# 300 digits); evaluations from several threads agree bit for bit with serial ones.
using Test, KerrGeodesics
const KG = KerrGeodesics

@testset "Float64 is unchanged by the precision helpers" begin
    for tol in (1.0e-6, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11, 2.0e-11, 1.0e-12, 2.0e-12, 1.0e-13,
            1.0e-14, 1.0e-15, 1.0e-18)
        @test KG._tol(Float64, tol) === tol
    end
    for n in (4, 10, 12, 32, 33, 64)
        @test KG._nterms(Float64, n) === n
    end
    T = Float64
    @test eps(T) === eps() && floatmin(T) === floatmin() && prevfloat(one(T)) === prevfloat(1.0)
    @test T(π) / 2 === pi / 2 && T(π) / 4 === π / 4 && 2T(π) === 2pi && T(π) === 1.0pi
    @test T(π) / 2 === 0.5pi && T(π) === acos(-1.0)
    @test sqrt(T(2)) === sqrt(2.0) && log(T(2)) === log(2) && 2log(T(2)) === 2log(2.0)
    @test T(2) / 3 === 2 / 3 && T(1) / 3 === 1 / 3 && inv(T(6)) === 1 / 6 && T(1) / 6 === 1 / 6
    @test T(3) / 40 === 3 / 40 && T(11) / 20 === 0.55 && T(2) / 5 === 2 / 5
    @test KG._float_type(1, 2.0) === Float64 && KG._float_type(big(1), 0.5) === BigFloat
    @test KG._real_type(1.0 + 2im, 3) === Float64
    @test KG._with_precision(() -> precision(BigFloat), Float64, 512) == precision(BigFloat)
    @test KG._with_precision(() -> precision(BigFloat), BigFloat, 512) == 512
end

@testset "Tolerances and truncation orders scale with the precision" begin
    setprecision(BigFloat, 256) do
        T = BigFloat
        for tol in (1.0e-8, 1.0e-12, 1.0e-14)
            t = KG._tol(T, tol)
            # the same place on the logarithmic scale between eps and 1
            @test log(t) / log(eps(T)) ≈ log(tol) / log(eps(Float64)) rtol = 1e-12
        end
        @test KG._nterms(T, 32) == cld(32 * 256, 53)
    end
end

# mpmath references: func, arguments (binary64, hex), value (170 digits), error scale
const ELLIPTIC_REFS = [split(l, '\t') for l in eachline(joinpath(@__DIR__, "elliptic_refs.tsv"))]
_jacobi(m) = KG._jacobi_parameter(m)
const ELLIPTIC_IMPL = Dict(
    "K" => (m,) -> KG._K(m), "E" => (m,) -> KG._E(m), "D" => (m,) -> KG._D(m),
    "Pi" => (n, m) -> KG._Pi(n, m), "F" => (p, m) -> KG._F(p, m), "Ephi" => (p, m) -> KG._E(p, m),
    "Dphi" => (p, m) -> KG._D(p, m), "Piphi" => (n, p, m) -> KG._Pi(n, p, m),
    "am" => (u, m) -> KG._am(u, _jacobi(m)), "sn" => (u, m) -> KG._sn(u, _jacobi(m)),
    "cn" => (u, m) -> KG._cn(u, _jacobi(m)), "dn" => (u, m) -> KG._dn(u, _jacobi(m)))

# largest error |value − reference|/scale over the references, at `bits` of precision
function elliptic_error(bits)
    worst = 0.0
    for (f, args, refstr, scalestr) in ELLIPTIC_REFS
        x = [parse(Float64, a) for a in split(args, ',')]
        ref, scale = setprecision(() -> (parse(BigFloat, refstr), parse(BigFloat, scalestr)), BigFloat, 700)
        v = bits == 53 ? big(ELLIPTIC_IMPL[f](x...)) :
            setprecision(() -> ELLIPTIC_IMPL[f](BigFloat.(x)...), BigFloat, bits)
        @test v isa BigFloat || bits == 53
        err = setprecision(BigFloat, 700) do
            isinf(ref) ? (v == ref ? 0.0 : Inf) :
                Float64(iszero(scale) ? abs(v) : abs(v - ref) / scale)
        end
        worst = max(worst, err)
    end
    return worst
end

@testset "Elliptic integrals and Jacobi functions in the precision of their arguments" begin
    @test elliptic_error(53) <= 16 * eps(Float64)
    @test elliptic_error(128) <= 16 * 2.0^-127
    @test elliptic_error(256) <= 1e-74
    @test elliptic_error(512) <= 1e-150
end

@testset "Elliptic domains" begin
    @test KG._K(1.0) == Inf && KG._E(1.0) == 1 && isnan(KG._K(NaN))
    @test_throws DomainError KG._K(1.5)
    @test_throws DomainError KG._F(0.3, -0.1)
    @test_throws DomainError KG._Pi(0.2, 1.1)
    @test_throws DomainError KG._am(0.4, KG._jacobi_parameter(1.2))
    @test KG._F(0.5, 1.0) ≈ atanh(sin(0.5)) && KG._F(-0.5, 1.0) ≈ -atanh(sin(0.5))
    @test all(isnan, KG._jacobi_sncndn(0.3, KG._jacobi_parameter(NaN)))
    @test KG._jacobi_sncndn(0.0, KG._jacobi_parameter(NaN)) == (0.0, 1.0, 1.0)
    @test KG._jacobi_sncndn(0.7, KG._jacobi_parameter(1.0)) == (tanh(0.7), sech(0.7), sech(0.7))
end

include("floats.jl")
include("roots.jl")
include("members.jl")
include("interfaces.jl")
include("axis.jl")
include("apex_geometry.jl")
include("threads.jl")
