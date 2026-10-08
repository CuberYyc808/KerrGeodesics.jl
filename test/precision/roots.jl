# Roots of the radial potential and the classification in the precision of the constants:
# every root of R (mpmath references at 300 digits), and the classification's numbers in T.

const ROOT_REFS = [split(l, '\t'; keepempty=true) for l in eachline(joinpath(@__DIR__, "root_refs.tsv"))]

# largest |root − reference|/max(1, |reference|) over all roots of all cases, at `bits`
function root_error(bits)
    worst = 0.0
    for (inputs, reals, pairs) in ROOT_REFS
        x = [parse(Float64, s) for s in split(inputs, ',')]
        refs = setprecision(BigFloat, 400) do
            [[complex(parse(BigFloat, r)) for r in split(reals, ',') if !isempty(r)];
             [complex(parse(BigFloat, split(p, ':')[1]), s * parse(BigFloat, split(p, ':')[2]))
              for p in split(pairs, ',') if !isempty(p) for s in (1, -1)]]
        end
        T = bits == 53 ? Float64 : BigFloat
        structure = setprecision(() -> kerr_geo_root_structure(T.(x)...), BigFloat, max(bits, 64))
        @test eltype(structure.raw_roots) == Complex{T}
        @test all(r -> r.radius isa T, structure.real_roots)
        found = structure.raw_roots
        @test length(found) == length(refs)
        for z in refs
            err = setprecision(BigFloat, 400) do
                Float64(minimum(abs(w - z) for w in found) / max(1, abs(z)))
            end
            worst = max(worst, err)
        end
    end
    return worst
end

@testset "Radial roots in the precision of the constants" begin
    @test root_error(53) <= 1e-14
    @test root_error(128) <= 1e-33
    @test root_error(256) <= 1e-70
end

@testset "Classification in BigFloat" begin
    setprecision(BigFloat, 256) do
        for (a, E, L, Q) in ((0.9, 0.94, 2.0, 4.0), (0.5, 1.1, 3.0, 2.0), (0.9, 0.97, 2.8, 5.0))
            c = kerr_geo_classify(big(a), big(E), big(L), big(Q))
            f = kerr_geo_classify(a, E, L, Q)
            @test c isa KerrGeoClassification{BigFloat}
            @test isempty(non_bigfloat(c))
            @test isempty(non_bigfloat(kerr_geo_root_structure(big(a), big(E), big(L), big(Q))))
            @test c.CaseIds == f.CaseIds
            @test all(k -> k.LowerEndpoint.Radius isa BigFloat, c.Components)
            @test all(k -> isapprox(k.UpperEndpoint.Radius, f.Components[1].UpperEndpoint.Radius;
                rtol=1e-12) || k !== c.Components[1], c.Components)
        end
    end
end

@testset "Chebyshev tables in the precision of their abscissae" begin
    for bits in (128, 256)
        setprecision(BigFloat, bits) do
            T = BigFloat
            f(x) = (exp(x), sin(3x), one(x) / (2 + x))
            p = KerrGeodesics.chebfit(f, [zero(T), T(1) / 2, T(2)]; ncomp=3)
            prim = KerrGeodesics.chebintegrate(p)
            @test p isa KerrGeodesics.ChebPieces{BigFloat}
            x = T(7) / 5
            # the tables meet the tolerance `_tol(T, 1e-14)` of the Float64 tables (1e-14)
            tol = 16 * KerrGeodesics._tol(T, 1e-14)
            @test abs(p(x, 1) - exp(x)) <= tol * exp(x)
            @test abs(prim(x, 1) - expm1(x)) <= tol * exp(x)
            @test abs(prim(x, 2) - (1 - cos(3x)) / 3) <= tol
            @test abs(prim(x, 3) - log((2 + x) / 2)) <= tol
            @test KerrGeodesics._tol(T, 1e-14) < (bits == 256 ? 1e-67 : 1e-31)
        end
    end
    p64 = KerrGeodesics.chebfit(x -> (exp(x),), 0.0, 1.0)
    @test length(p64.coefs[1][1]) <= KerrGeodesics._cheb_n(Float64) + 1
    @test KerrGeodesics._cheb_n(Float64) == 32
end

@testset "Far-field r* and horizon azimuth in the working precision" begin
    setprecision(BigFloat, 256) do
        # u = r − 1 > 200 (r₊ − r₋): the expansion in r₊ − r₋ against the closed form in 512 bits
        for (a, r) in ((big"0.6", big(400)), (big"0.999", big(30)), (big"0.3", big(1000)))
            exact = setprecision(BigFloat, 512) do
                A, R = BigFloat(a, precision=512), BigFloat(r, precision=512)
                rp, rm = 1 + sqrt(1 - A^2), 1 - sqrt(1 - A^2)
                d = rp - rm
                (R + 2rp / d * log((R - rp) / 2) - 2rm / d * log((R - rm) / 2),
                 A / d * log((R - rp) / (R - rm)))
            end
            @test abs(kerr_rstar(a, r) - exact[1]) <= 1e-70 * abs(exact[1])
            @test abs(KerrGeodesics._horizon_azimuth(a, r) - exact[2]) <= 1e-70 * abs(exact[2])
        end
    end
end
