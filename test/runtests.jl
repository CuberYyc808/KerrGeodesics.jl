using Test
using KerrGeodesics

@testset "basic" begin
    # a tiny smoke test;
    vals = Kerr_Geodesics(0.0, 10.0, 0.5, 0.8)
    @test isa(vals, Dict) || isa(vals, NamedTuple)
end
