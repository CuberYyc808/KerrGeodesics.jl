using Test
using KerrGeodesics

@testset "KerrGeodesics.jl" begin
    plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:inner_turning)
    @test plunge.Status.supported
    @test plunge.OrbitClass == "Complex"
    @test isapprox(plunge.Trajectory.r(0.0), 0.4963292130784854; atol=1e-12)
    @test abs(plunge.Residuals.radial(0.5)) < 1e-10
    @test abs(plunge.Residuals.polar_z(0.5)) < 1e-10

    exterior = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)
    @test isfinite(exterior.Trajectory.rstar(0.0))
    @test isapprox(exterior.Trajectory.u(0.0), exterior.Trajectory.t(0.0) - exterior.Trajectory.rstar(0.0); atol=1e-12)
    @test isapprox(exterior.Trajectory.v(0.0), exterior.Trajectory.t(0.0) + exterior.Trajectory.rstar(0.0); atol=1e-12)
    @test isnan(plunge.Trajectory.rstar(0.0))

    direct = kerr_geodesic(0.9, (0.94, 0.1, 12.0); radial_start=:turning_point)
    @test direct.InputType == :constants
    @test direct isa KerrGeodesicFamily
    @test direct isa KerrGeodesicS
    @test direct isa KerrGeodesicSet
    @test direct.Plunge.Status.supported
    @test direct.Plunge.Status.duration.mino_time_to_horizon > 0
    @test direct.Plunge.Status.duration.coordinate_time_to_horizon_status == "boyer_lindquist_t_diverges_at_future_horizon"
    @test direct.Stable === nothing

    combined = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
    @test combined.InputType == :apex
    @test combined.Stable !== nothing
    @test combined.Plunge.Status.supported

    near_sep_p = kerr_geo_orbit_type_metadata(0.9, 10.0, 0.5, 0.8).separatrix_p + 1e-11
    near_sep_type = kerr_geo_orbit_type_metadata(0.9, near_sep_p, 0.5, 0.8)
    @test near_sep_type.family == "Bound"
    @test near_sep_type.start_type == "NotApplicable"
    @test near_sep_type.stability in ("Stable", "MarginallyStable")
    @test !("BoundPlunge" in near_sep_type.labels)
    near_sep_freqs = kerr_geo_frequencies(0.9, near_sep_p, 0.5, 0.8; Time="Mino")
    @test all(isfinite, values(near_sep_freqs))
    @test near_sep_freqs["ϒr"] >= 0

    exact_sep_p = near_sep_type.separatrix_p
    sep_minus_roundoff = kerr_geo_orbit_type_metadata(0.9, exact_sep_p - 5e-16, 0.5, 0.8)
    sep_plus_roundoff = kerr_geo_orbit_type_metadata(0.9, exact_sep_p + 5e-16, 0.5, 0.8)
    @test sep_minus_roundoff.at_separatrix
    @test sep_plus_roundoff.at_separatrix
    @test sep_minus_roundoff.effective_p == sep_minus_roundoff.separatrix_p
    @test sep_plus_roundoff.effective_p == sep_plus_roundoff.separatrix_p
    @test "Separatrix" in sep_minus_roundoff.labels
    @test "Separatrix" in sep_plus_roundoff.labels
    @test !("MarginallyStable" in sep_minus_roundoff.labels)
    @test !("MarginallyStable" in sep_plus_roundoff.labels)
    @test !("Unstable" in sep_minus_roundoff.labels)
    @test sep_minus_roundoff.family == "Bound"
    @test !("BoundPlunge" in sep_minus_roundoff.labels)

    bound_plunge_type = kerr_geo_orbit_type_metadata(0.9, near_sep_type.separatrix_p - 1e-5, 0.5, 0.8)
    @test bound_plunge_type.family == "Plunge"
    @test bound_plunge_type.start_type == "BoundPlunge"
    @test bound_plunge_type.support_status == "bound_plunge_metadata_only"
    @test !("Unstable" in bound_plunge_type.labels)

    scatter_type = kerr_geo_orbit_type_metadata(0.9, kerr_geo_separatrix(0.9, 1.2, 0.8) + 1e-5, 1.2, 0.8)
    @test scatter_type.family == "Scatter"
    @test scatter_type.start_type == "InfinityStart"
    @test scatter_type.support_status == "scattering_metadata_only_not_implemented"

    infinity_plunge_type = kerr_geo_orbit_type_metadata(0.9, kerr_geo_separatrix(0.9, 1.2, 0.8) - 1e-5, 1.2, 0.8)
    @test infinity_plunge_type.family == "Plunge"
    @test infinity_plunge_type.start_type == "InfinityStart"
    @test infinity_plunge_type.support_status == "infinity_start_plunge_metadata_only_not_implemented"

    real2 = kerr_geo_plunge(
        0.86, 0.66, 0.85, 1e-2;
        initial_radius=2.8722246238500517,
        initial_theta=pi / 2,
    )
    @test real2.Status.supported
    @test real2.OrbitClass == "Real2"
    @test real2.Status.time_phi_status == "stable_realroot_analytic_branch_continuation_with_linear_anchor"
    @test real2.Status.rstar_convention == "code_log_halves"
    @test real2.Status.near_horizon_branch == "ingoing_future_horizon"
    @test real2.Status.P_plus_sign == "positive_required"
    @test real2.Status.near_horizon_series_order == 10
    @test real2.Status.u_status == "diverges_linearly_for_ingoing_future_horizon"
    @test real2.Status.v_status == "finite_regular_horizon_coordinate"
    @test real2.Status.lambda_to_rstar_path == "numerical_or_existing_trajectory_not_replaced"
    rstar_end = real2.Trajectory.rstar(real2.Status.real2_lambda_end)
    @test isfinite(real2.Trajectory.v_rstar_series(rstar_end))
    @test isapprox(
        real2.Trajectory.u_rstar_series(rstar_end),
        real2.Trajectory.v_rstar_series(rstar_end) - 2 * rstar_end;
        atol=1e-12,
    )
    @test isfinite(real2.Trajectory.t(0.1))
    @test isfinite(real2.Trajectory.phi(0.1))
    @test abs(real2.Residuals.radial(0.1)) < 1e-10
    @test abs(real2.Residuals.polar_z(0.1)) < 1e-10
end
