# KerrGeodesics.jl

Julia interfaces for Kerr geodesic trajectories in units with `G = c = M = 1`.

Current support:

- stable bound orbits from APEX-like parameters `(a,p,e,x)`;
- bound plunge orbits from constants `(a,E,Lz,Q)`, with `E < 1`.

Scattering orbits are planned for future development.

## Installation

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

## Basic Usage

```julia
using KerrGeodesics

stable = kerr_geo_stable(0.9, 10.0, 0.5, 0.8)
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)
family = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
```

Typical printed outputs are:

```julia
KerrGeoStable(
    OrbitalParameters = (a = 0.9, p = 10.0, e = 0.5, x = 0.8),
    ConstantsOfMotion = (E = 0.9641204328952226, Lz = 2.8359152778998453, Q = 4.544408272395823),
    OrbitalType = ["Bound", "Eccentric", "Stable", "Inclined"],
    Frequencies = (ϒt = 171.0926187383033, ϒr = 2.792721794117058, ϒθ = 3.551489601048812, ϒϕ = 3.7357605214030265),
    Parametrization = "Mino",
    Trajectory = (t = t(λ), r = r(λ), θ = θ(λ), ϕ = ϕ(λ)),
    InitialPhases = (qt0 = 0.0, qr0 = 0.0, qθ0 = 0.0, qϕ0 = 0.0),
)
```

```julia
KerrGeoPlunge(
    ConstantsOfMotion = (E = 0.94, Lz = 0.1, Q = 12.0),
    OrbitClass = "Complex",
    Parametrization = "Mino",
    InitialPosition = (t0 = 0.0, r0 = 3.203955290691315, theta0 = 1.5707963267948966, phi0 = 0.0),
    Status = (supported = true, reason = "ok", trajectory_support_status = "production_trajectory_supported", production_branch_support = "supported_complex_current_route", time_phi_status = "analytic_branch_formula", mathematica_reference_status = "available_complex_reference_record", source_route_status = "supported_complex_current_route", waveform_handoff_status = "metadata_only_unless_source_route_validated", retarded_advanced_time_status = "trajectory_callables_available_finite_cutoff_convention_required_near_horizon", real2_horizon_offset = NaN, real2_rplus = NaN, real2_r_end = NaN, real2_lambda_end = NaN, rstar_convention = "code_log_halves", near_horizon_branch = "unavailable_for_current_branch", P_plus_sign = "not_evaluated", near_horizon_series_order = 0, near_horizon_series_variable = "q_exp_rstar_minus_Cstar_over_alpha", u_status = "direct_trajectory_callable_only", v_status = "direct_trajectory_callable_only", lambda_to_rstar_path = "numerical_or_existing_trajectory_not_replaced", series_switch_tolerance_target = 1.0e-15, series_switch_status = "not_configured", regular_time_coefficients_status = "not_configured", time_origin = :input_t0, horizon_time_origin_status = "input_t0_preserved", horizon_time_shift = 0.0, horizon_v_anchor_before_shift = NaN, horizon_v_anchor_method = "input_t0_no_horizon_alignment", horizon_time_origin_offset = 0.0001, duration = (start_lambda = 0.0, horizon_lambda = 1.116067929320647, mino_time_to_horizon = 1.116067929320647, radial_endpoint_lambda = 1.6863825015685567, mino_time_to_radial_endpoint = 1.6863825015685567, horizon_radius = 1.4358898943540672, coordinate_time_to_horizon_status = "boyer_lindquist_t_diverges_at_future_horizon", retarded_time_to_horizon_status = "u_diverges_linearly_at_future_horizon", advanced_time_to_horizon_status = "v_has_finite_horizon_anchor"), horizon_lambda = 1.116067929320647, mino_time_to_horizon = 1.116067929320647, radial_endpoint_lambda = 1.6863825015685567, mino_time_to_radial_endpoint = 1.6863825015685567),
    Trajectory = (t = t(lambda), r = r(lambda), theta = theta(lambda), phi = phi(lambda), rstar = rstar(lambda), u = u(lambda), v = v(lambda), u_rstar_series = u(rstar), v_rstar_series = v(rstar)),
    Velocity = (ut = ut(lambda), ur = ur(lambda), uz = dz/dlambda, utheta = dtheta/dlambda, uphi = uphi(lambda)),
)
```

```julia
KerrGeodesicFamily(
    InputType = :apex,
    Parameters = (a = 0.9, p = 10.0, e = 0.5, x = 0.8),
    ConstantsOfMotion = (E = 0.9641204328952226, Lz = 2.8359152778998453, Q = 4.544408272395823),
    RootClass = "Real1",
    HasStable = true,
    HasPlunge = true,
    Status = (supported = true, reason = "ok"),
)
```

The family object stores the two compatible orbit objects directly:

```julia
family.Stable
family.Plunge
```

For example:

```julia
family.Stable.Trajectory.r(0.0)
# 6.666666666666667

family.Plunge.OrbitClass
# "Real1"

family.Plunge.Status.duration.mino_time_to_horizon
# 0.05426602766253312
```

For a stable eccentric orbit, `initPhases=(0,0,0,0)` starts the radial motion at
periapsis, `r(0)=p/(1+e)`. For a bound plunge, `radial_start=:turning_point`
starts at the exterior turning point.

## Examples

The example notebook builds both stable and bound-plunge trajectory animations:

- [`example/Test_KerrGeodesics.ipynb`](example/Test_KerrGeodesics.ipynb)
- [`example/generate_example_gifs.jl`](example/generate_example_gifs.jl)

Generated example images:

![Stable bound Kerr geodesic](example/Trajectory_stable.gif)

![Bound plunge Kerr geodesic](example/Trajectory_plunge.gif)

## Documentation

Documenter.jl sources live under [`docs/`](docs/).

## License

MIT.
