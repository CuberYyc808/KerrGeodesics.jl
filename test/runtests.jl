# Published test: the package loads and answers one call. The full test suite is kept in the
# development repository.
using Test
using KerrGeodesics

@testset "KerrGeodesics loads" begin
    @test isdefined(KerrGeodesics, :kerr_geodesic)
    @test 0 < kerr_geo_constants_of_motion(0.9, 10.0, 0.5, 0.8)["E"] < 1
    family = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
    @test family.Stable isa KerrGeoStableComponent
    @test 0 < family.Stable.ConstantsOfMotion.E < 1
end

@testset "Compact default display" begin
    family = kerr_geodesic(0.9, (0.94, 0.1, 12.0))
    member = family.Plunge
    for object in (family, member)
        text = repr("text/plain", object)
        @test length(split(text, '\n')) <= 8
        @test occursin("Status", text) && occursin("supported", text)
        @test !occursin("spectral=", text) && !occursin("member_errors=", text)
        @test !occursin("ReferenceZero", text) && !occursin("#", text)
    end
    @test occursin("(t(lambda), r(lambda), theta(lambda), phi(lambda))",
        repr("text/plain", member))
    @test member.Status isa NamedTuple && member.Trajectory isa NamedTuple
    @test kerr_geo_members(family)[end] === member
end

include(joinpath(@__DIR__, "precision", "runtests.jl"))
