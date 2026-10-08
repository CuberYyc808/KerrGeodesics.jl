# APEX input keeps its turning points. Near the separatrix the rounded constants alone can merge
# r2 and r3 into a repeated root (classified K3/K4/K5, no Stable member); the Stable member is
# built from r1 = p/(1 − e), r2 = p/(1 + e) and the inner quadratic of the same constants, and
# every other member keeps the roots of the constants.

@testset "APEX turning points near the separatrix" begin
    # a = 0, e = 1/10, x = −1: separatrix p = 6 + 2e = 31/5; p = 31/5 + 1e-7 (exact)
    a, e, x = 0.0, 0.1, -1.0
    p = Float64(31 // 5 + 1 // 10^7)
    com = kerr_geo_constants_of_motion(a, p, e, x)
    constants = kerr_geodesic(a, (com["E"], com["Lz"], com["Q"]))
    @test constants.Status.classification.CaseIds == (:K3, :K4, :K5)
    @test constants.Stable === nothing
    @test !haskey(constants.Status, :apex_root_geometry)
    kg = kerr_geodesic(a, p, e, x)
    g = kg.Status.apex_root_geometry
    @test g.accepted && g.reason === :accepted
    @test kg.Stable !== nothing && kg.Stable.CaseId === :A1
    @test isempty(kg.Status.member_errors)
    @test kg.Stable.Roots.radial[1:2] == (p / (1 - e), p / (1 + e))
    # Schwarzschild: r3 = 2p/(p − 4), r4 = 0
    @test abs(kg.Stable.Roots.radial[3] - 2p / (p - 4)) <= 128eps(p)
    @test iszero(kg.Stable.Roots.radial[4])
    @test g.diagnostics.gap > g.diagnostics.gap_bound > 0
    @test max(g.diagnostics.root_backward, g.diagnostics.coefficient_backward) <= 256eps()
    # one root model per component, stated in the status
    @test kg.Status.case_ids == (:A1, :K3, :K4, :K5)
    @test all(m -> m.roots === (m.case_id === :A1 ? :apex_turning_points : :constants),
        kg.Status.component_root_models)
    @test Tuple(c.CaseId for c in kg.Critical) == (:K3, :K4, :K5)
    @test kerr_geodesic(a, p, e, x; case_id=:A1).Status.selected_case === :A1
    @test kerr_geodesic(a, p, e, x; case_id=:K3).Status.selected_case === :K3
    # the same orbit at 128 bits: the constants resolve the gap themselves
    kb = kerr_geodesic(0 // 1, 31 // 5 + 1 // 10^7, 1 // 10, -1 // 1; precision=128)
    @test kb.Status.classification.CaseIds == (:A1, :B1)
    @test kb.Status.apex_root_geometry.accepted
    @test abs(kb.Stable.Status.frequencies.ϒr - big(kg.Stable.Status.frequencies.ϒr)) <=
        1e-9 * kb.Stable.Status.frequencies.ϒr
end

@testset "APEX turning points: unresolved geometry keeps the constants" begin
    # near-parabolic e = 0.999: rejected by the backward-error contract
    kg = kerr_geodesic(0.5, 20.0, 0.999, 0.8)
    @test !kg.Status.apex_root_geometry.accepted
    @test kg.Status.apex_root_geometry.reason === :backward_contract
    @test all(m -> m.roots === :constants, kg.Status.component_root_models)
    # circular orbits (e = 0) and the separatrix itself are outside the geometry's domain
    @test kerr_geodesic(0.9, 10.0, 0.0, 0.8).Status.apex_root_geometry.reason === :outside_domain
    sep = kerr_geodesic(0.0, 7.0, 0.5, 1.0)
    @test !sep.Status.apex_root_geometry.accepted
end
