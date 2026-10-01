# Published test: the package loads and answers one call. The full test suite is kept in the
# development repository.
using Test
using KerrGeodesics

@testset "KerrGeodesics loads" begin
    @test isdefined(KerrGeodesics, :kerr_geodesic)
    @test 0 < kerr_geo_constants_of_motion(0.9, 10.0, 0.5, 0.8)["E"] < 1
end
