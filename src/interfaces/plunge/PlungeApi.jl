# Plunge reference API: the `KerrGeoPlunge` record and `kerr_geo_plunge`.

"""
    KerrGeoPlunge

Result of `kerr_geo_plunge`: an E < 1 Kerr plunge built from the spin and constants of
motion (a, E, Lz, Q). `OrbitClass` is the radial root class (Real1, Real2 or Complex);
`Trajectory`, `Velocity` and `Residuals` hold functions of Mino time λ, `Potentials` the
radial and polar potentials R(r) and Θ(z), and `Status` the conventions and the Mino time
to the horizon.
"""
struct KerrGeoPlunge
    OrbitClass::String
    OrbitalParameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    Parametrization::String
    Roots::NamedTuple
    InitialPhases::NamedTuple
    Trajectory::Any
    Velocity::Any
    Potentials::NamedTuple
    Residuals::NamedTuple
    Status::NamedTuple
end

function Base.show(io::IO, kg::KerrGeoPlunge)
    print(io, "KerrGeoPlunge(constants=")
    show(io, kg.ConstantsOfMotion)
    print(io, ", supported=", kg.Status.supported, ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoPlunge)
    println(io, "KerrGeoPlunge")
    _show_summary_field(io, "Parameters", kg.OrbitalParameters)
    _show_summary_field(io, "Constants", kg.ConstantsOfMotion)
    _show_summary_field(io, "Root class", kg.OrbitClass)
    println(io, "  Trajectory = (t(lambda), r(lambda), theta(lambda), phi(lambda))")
    _show_summary_status(io, kg.Status)
end

"""
    kerr_geo_plunge(a, E, Lz, Q; kwargs...)

Build an E < 1 plunge from the Kerr spin and constants of motion. `Trajectory` holds `t`,
`r`, `theta`, `phi`, `rstar`, `u = t − r*` and `v = t + r*` as functions of Mino time λ for
every root class (Real1, Real2, Complex). `OrbitClass` is the root class; `Status` records
the Mino time to the horizon (`duration`), the near-horizon series and the time origin.

Keywords: `radial_start = :turning_point` (default; λ = 0 at the plunge turning point
outside the horizon) or `:inner_turning` (λ = 0 at the inner end of the radial range);
`initial_radius` and `initial_theta` (default π/2) start the orbit at that radius on the
ingoing leg and at that polar angle; `radial_phase`, `theta_phase` or
`initPhases = (t0, λr0, λθ0, φ0)` give the phases directly; `t0` and `phi0` offset t and φ.
`time_origin = :input_t0` (default) keeps t(0) = t0; `:future_horizon_v_zero` shifts t so
that v = t + r* vanishes on the future horizon, with the horizon value of v taken from radii
`horizon_time_origin_offset` above r₊.

For the Real2 root class, `radial_start=:inner_turning` puts the start event (λ = 0) on the
future horizon r₊ itself. Boyer-Lindquist t and φ diverge logarithmically there, so t and φ
measured from that event are not defined: `t` and `phi` (and u = t − r*, v = t + r* built
from t) return NaN for every λ; r and θ remain available.
"""
function kerr_geo_plunge(a::Real, energy::Real, lz::Real, q::Real;
        initPhases=nothing,
        radial_phase=nothing,
        theta_phase=nothing,
        initial_radius=nothing,
        initial_theta=pi / 2,
        radial_start=:turning_point,
        t0=0.0,
        phi0=0.0,
        real2_horizon_offset=1e-4,
        time_origin=:input_t0,
        horizon_time_origin_offset=real2_horizon_offset)

    roots, root_class = classify_orbit(a, energy, lz, q)

    zm, zp = polar_roots(a, energy, lz, q)
    if initPhases === nothing
        lambda_r0 = _radial_phase_from_options(
            a, energy, lz, q;
            radial_phase=radial_phase,
            initial_radius=initial_radius,
            radial_start=radial_start,
        )
        lambda_theta0 = theta_phase === nothing ?
            _theta_phase_from_theta(a, energy, lz, q, initial_theta) :
            theta_phase
        phases = (t0, lambda_r0, lambda_theta0, phi0)
    else
        phases = initPhases
        lambda_r0 = phases[2]
        lambda_theta0 = phases[3]
    end

    t_bl, r, theta, phi_bl = generic_plunge_orbit(
        a, energy, lz, q;
        initPhases=phases,
        real2_horizon_offset=real2_horizon_offset,
    )
    lambda_radial_endpoint, lambda_horizon_from_turning_point, lambda_of_radius = lambda_of_r(a, energy, lz, q)
    rplus = _rplus(a)
    horizon_lambda_from_start = lambda_horizon_from_turning_point - lambda_r0
    # a start event on the horizon (Real2 :inner_turning) has no Boyer-Lindquist t, φ
    on_horizon = iszero(horizon_lambda_from_start)
    t_raw = on_horizon ? (λ -> NaN) : t_bl
    phi = on_horizon ? (λ -> NaN) : phi_bl
    radial_endpoint_lambda_from_start = lambda_radial_endpoint - lambda_r0
    duration_metadata = (
        start_lambda=0.0,
        horizon_lambda=horizon_lambda_from_start,
        mino_time_to_horizon=horizon_lambda_from_start,
        radial_endpoint_lambda=radial_endpoint_lambda_from_start,
        mino_time_to_radial_endpoint=radial_endpoint_lambda_from_start,
        horizon_radius=rplus,
    )
    time_origin_anchor = (
        vH=NaN,
        method="input_t0_no_horizon_alignment",
        direct_anchor_values=Float64[],
    )
    horizon_time_shift = 0.0
    if time_origin == :input_t0
        horizon_time_shift = 0.0
    elseif time_origin == :future_horizon_v_zero
        time_origin_anchor = _estimate_future_horizon_v_anchor(
            a, energy, lz, q, root_class, theta, t_raw, lambda_of_radius, lambda_r0;
            order=10,
            horizon_offset=horizon_time_origin_offset,
        )
        horizon_time_shift = -time_origin_anchor.vH
    else
        error("Unknown time_origin. Use :input_t0 or :future_horizon_v_zero.")
    end
    t(lambda) = t_raw(lambda) + horizon_time_shift
    phases_effective = (
        phases[1] + horizon_time_shift,
        phases[2],
        phases[3],
        phases[4],
    )
    ut, ur, uz, uphi = generic_plunge_velocity(a, energy, lz, q; initPhase=(lambda_r0, lambda_theta0))
    utheta(lambda) = begin
        s = sin(theta(lambda))
        abs(s) < sqrt(eps(Float64)) ? NaN : -uz(lambda) / s
    end
    rstar(lambda) = kerr_rstar(a, r(lambda))
    u(lambda) = t(lambda) - rstar(lambda)
    v(lambda) = t(lambda) + rstar(lambda)

    radial_potential(rvalue) = kerr_radial_potential(a, energy, lz, q, rvalue)
    polar_potential_z(zvalue) = kerr_polar_z_potential(a, energy, lz, q, zvalue)
    radial_residual(lambda) = ur(lambda)^2 - radial_potential(r(lambda))
    polar_residual(lambda) = uz(lambda)^2 - polar_potential_z(cos(theta(lambda)))
    real2_rplus = root_class == "Real2" ? 1 + sqrt(1 - a^2) : NaN
    real2_r_end = root_class == "Real2" ? real2_rplus + real2_horizon_offset : NaN
    real2_lambda_end = if root_class == "Real2"
        lambda_of_radius(real2_r_end) - lambda_r0
    else
        NaN
    end
    near_horizon_series = if root_class == "Real2"
        _build_near_horizon_uv_series(
            a, energy, lz, q, theta, t, lambda_of_radius, lambda_r0;
            order=10,
            horizon_offset=real2_horizon_offset,
        )
    else
        _near_horizon_series_unavailable()
    end
    near_horizon_metadata = near_horizon_series.metadata

    real2 = root_class == "Real2"
    near_horizon_fields = NamedTuple{(:rstar_convention, :near_horizon_branch, :P_plus_sign,
        :near_horizon_series_order, :near_horizon_series_variable,
        :lambda_to_rstar_path)}(near_horizon_metadata)
    status = (; supported=true, reason="ok",
        real2_horizon_offset=real2 ? real2_horizon_offset : NaN,
        real2_rplus=real2_rplus, real2_r_end=real2_r_end, real2_lambda_end=real2_lambda_end,
        near_horizon_fields...,
        time_origin=time_origin,
        horizon_time_shift=horizon_time_shift,
        horizon_v_anchor_before_shift=time_origin_anchor.vH,
        horizon_v_anchor_method=time_origin_anchor.method,
        horizon_time_origin_offset=horizon_time_origin_offset,
        duration=duration_metadata,
        horizon_lambda=horizon_lambda_from_start,
        mino_time_to_horizon=horizon_lambda_from_start,
        radial_endpoint_lambda=radial_endpoint_lambda_from_start,
        mino_time_to_radial_endpoint=radial_endpoint_lambda_from_start)

    return KerrGeoPlunge(
        root_class,
        (a=a,),
        (E=energy, Lz=lz, Q=q),
        "Mino",
        (radial=roots, polar=(zm=zm, zp=zp)),
        (t0=phases_effective[1], radial=phases_effective[2], theta=phases_effective[3], phi0=phases_effective[4]),
        (
            t=t,
            r=r,
            theta=theta,
            phi=phi,
            rstar=rstar,
            u=u,
            v=v,
            u_rstar_series=near_horizon_series.u_rstar,
            v_rstar_series=near_horizon_series.v_rstar,
            near_horizon_q=near_horizon_series.q_of_rstar,
            near_horizon_series_last_term_abs=near_horizon_series.last_term_abs,
        ),
        (ut=ut, ur=ur, uz=uz, utheta=utheta, uphi=uphi),
        (radial=radial_potential, polar_z=polar_potential_z),
        (radial=radial_residual, polar_z=polar_residual),
        status
    )
end
