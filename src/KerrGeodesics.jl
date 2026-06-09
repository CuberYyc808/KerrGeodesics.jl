module KerrGeodesics

include("KerrGeoStable/ConstantsOfMotion.jl")
using .ConstantsOfMotion
include("KerrGeoStable/FourVelocity.jl")
using .FourVelocity
include("KerrGeoStable/OrbitalFrequencies.jl")
using .OrbitalFrequencies
include("KerrGeoStable/SpecialOrbits.jl")
using .SpecialOrbits
include("KerrGeoStable/KerrGeoOrbit.jl")
using .KerrGeoOrbit

export kerr_geo_constants_of_motion,
        kerr_geo_four_velocity,
        kerr_geo_frequencies,
        kerr_geo_radial_roots,
        kerr_geo_polar_roots,
        kerr_geo_orbit_type,
        kerr_geo_orbit_type_metadata,
        kerr_geo_separatrix,
        kerr_geo_isco,
        kerr_geo_ibso,
        kerr_geo_isso,
        kerr_geo_orbit,
        kerr_geo_stable,
        KerrGeoStable

"""
    KerrGeoStable

Structured output for a stable bound Kerr geodesic constructed from APEX-like
parameters `(a,p,e,x)`. The trajectory and four-velocity fields are callable
functions of Mino time.
"""
struct KerrGeoStable
    OrbitalType::Vector{String}
    OrbitalParameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    Parametrization::String
    Trajectory::NamedTuple
    InitialPhases::NamedTuple
    FourVelocity::NamedTuple
    Frequencies::NamedTuple
    CrossFunctions::NamedTuple
    DCrossFunctions::NamedTuple
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoStable)
    println(io, "KerrGeoStable(")
    print(io, "    OrbitalParameters = "); show(io, kg.OrbitalParameters); println(io, ",")
    print(io, "    ConstantsOfMotion = "); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    OrbitalType = "); show(io, kg.OrbitalType); println(io, ",")
    print(io, "    Frequencies = "); show(io, kg.Frequencies); println(io, ",")
    print(io, "    Parametrization = "); show(io, kg.Parametrization); println(io, ",")
    print(io, "    Trajectory = (t = t(λ), r = r(λ), θ = θ(λ), ϕ = ϕ(λ))"); println(io, ",")
    print(io, "    InitialPhases = "); show(io, kg.InitialPhases); println(io, ",")
    print(io, ")")
end

"""
    kerr_geo_stable(a, p, e, x; initPhases=(0.0, 0.0, 0.0, 0.0))

Construct a stable bound Kerr geodesic from APEX-like parameters. The phase
tuple is `(qt0, qr0, qtheta0, qphi0)`.
"""
function kerr_geo_stable(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    # Orbital Type
    otype = kerr_geo_orbit_type(a, p, e, x)

    # Constants of Motion
    com = kerr_geo_constants_of_motion(a, p, e, x)
    En = com["E"]
    L = com["Lz"]
    Q = com["Q"]

    # Trajectory
    KG = kerr_geo_orbit(a, p, e, x; initPhases = initPhases)
    t, r, θ, ϕ = KG["Trajectory"]
    # Frequencies
    freqs = kerr_geo_frequencies(a, p, e, x; Time="Mino")
    ϒt = freqs["ϒt"]
    ϒr = freqs["ϒr"]
    ϒθ = freqs["ϒθ"]
    ϒϕ = freqs["ϒϕ"]
    # Cross functions
    if KG["CrossFunction"] !== nothing
        Δtr = KG["CrossFunction"][1]
        Δtθ = KG["CrossFunction"][2]
        Δϕr = KG["CrossFunction"][3]
        Δϕθ = KG["CrossFunction"][4]
    else
        Δtr = nothing
        Δtθ = nothing
        Δϕr = nothing
        Δϕθ = nothing
    end
    # Derivatives of cross functions
    if KG["DerivativesCrossFunction"] !== nothing
        dtr = KG["DerivativesCrossFunction"][1]
        dtθ = KG["DerivativesCrossFunction"][2]
        dϕr = KG["DerivativesCrossFunction"][3]
        dϕθ = KG["DerivativesCrossFunction"][4]
    else
        dtr = nothing
        dtθ = nothing
        dϕr = nothing
        dϕθ = nothing
    end
    # Four-velocity
    ut, ur, uθ, uϕ = KG["FourVelocity"]
    
    return KerrGeoStable(
        otype,
        (a=a, p=p, e=e, x=x),
        (E=En, Lz=L, Q=Q),
        "Mino",
        (t=t, r=r, θ=θ, ϕ=ϕ),
        (qt0 = initPhases[1], qr0=initPhases[2], qθ0=initPhases[3], qϕ0=initPhases[4]),
        (ut=ut, ur=ur, uθ=uθ, uϕ=uϕ),
        (ϒt=ϒt, ϒr=ϒr, ϒθ=ϒθ, ϒϕ=ϒϕ),
        (Δtr=Δtr, Δtθ=Δtθ, Δϕr=Δϕr, Δϕθ=Δϕθ),
        (dtr=dtr, dtθ=dtθ, dϕr=dϕr, dϕθ=dϕθ)
    )
end


include("KerrGeoPlunge/OrbitClass.jl")
using .OrbitClass
include("KerrGeoPlunge/OrbitalDuration.jl")
using .OrbitalDuration
include("KerrGeoPlunge/FourVelocity.jl")
using .PlungeFourVelocity
include("KerrGeoPlunge/PlungeOrbit.jl")
using .PlungeOrbit
include("KerrGeoPlunge/NearHorizonTime.jl")
using .NearHorizonTime
include("KerrGeoPlunge/InitialConditions.jl")
using .InitialConditions

export radial_roots, polar_roots, classify_orbit, 
        lambda_of_r, 
        generic_plunge_velocity,
        generic_plunge_orbit,
        kerr_geo_orbit_type_metadata,
        kerr_rstar,
        KerrGeoPlunge,
        KerrGeodesicFamily,
        KerrGeodesicS,
        KerrGeodesicSet,
        kerr_geo_plunge,
        kerr_geodesic

"""
    KerrGeoPlunge

Structured output for a bound plunge Kerr plunge trajectory constructed from
constants of motion `(a,E,Lz,Q)`. Supported branches expose callable trajectory,
velocity, potential, residual, and metadata fields.
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

"""
    KerrGeodesicFamily

Combined project-facing output from `kerr_geodesic`. A family can contain a
stable branch, a plunge branch, or both, depending on the supplied parameter
type and root classification. `KerrGeodesicS` and `KerrGeodesicSet` are retained
as compatibility aliases.
"""
struct KerrGeodesicFamily
    InputType::Symbol
    Parameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    RootClass::String
    Stable::Any
    Plunge::Any
    Status::NamedTuple
end

const KerrGeodesicS = KerrGeodesicFamily
const KerrGeodesicSet = KerrGeodesicFamily

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoPlunge)
    println(io, "KerrGeoPlunge(")
    print(io, "    ConstantsOfMotion = "); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    OrbitClass = "); show(io, kg.OrbitClass); println(io, ",")
    print(io, "    Parametrization = "); show(io, kg.Parametrization); println(io, ",")
    initial_position = kg.Status.supported ?
        (t0=kg.InitialPhases.t0, r0=kg.Trajectory.r(0.0), theta0=kg.Trajectory.theta(0.0), phi0=kg.InitialPhases.phi0) :
        (t0=kg.InitialPhases.t0, r0=NaN, theta0=NaN, phi0=kg.InitialPhases.phi0)
    print(io, "    InitialPosition = "); show(io, initial_position); println(io, ",")
    if kg.Status.supported
        println(io, "    Trajectory = (t = t(lambda), r = r(lambda), theta = theta(lambda), phi = phi(lambda), rstar = rstar(lambda), u = u(lambda), v = v(lambda), u_rstar_series = u(rstar), v_rstar_series = v(rstar)),")
        println(io, "    Velocity = (ut = ut(lambda), ur = ur(lambda), uz = dz/dlambda, utheta = dtheta/dlambda, uphi = uphi(lambda)),")
    else
        println(io, "    Trajectory = nothing,")
        println(io, "    Velocity = nothing,")
    end
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeodesicFamily)
    println(io, "KerrGeodesicFamily(")
    print(io, "    InputType = "); show(io, kg.InputType); println(io, ",")
    print(io, "    Parameters = "); show(io, kg.Parameters); println(io, ",")
    print(io, "    ConstantsOfMotion = "); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    RootClass = "); show(io, kg.RootClass); println(io, ",")
    print(io, "    HasStable = "); show(io, kg.Stable !== nothing); println(io, ",")
    print(io, "    HasPlunge = "); show(io, kg.Plunge !== nothing && kg.Plunge.Status.supported); println(io, ",")
    print(io, "    Status = "); show(io, kg.Status); println(io, ",")
    print(io, ")")
end

function _radial_potential(a, energy, lz, q, r)
    return (energy * (r^2 + a^2) - a * lz)^2 -
           (r^2 - 2 * r + a^2) * (r^2 + (a * energy - lz)^2 + q)
end

function _polar_potential_z(a, energy, lz, q, z)
    return q * (1 - z^2) -
           z^2 * (a^2 * (1 - energy^2) * (1 - z^2) + lz^2)
end

function _unsupported_plunge(a, energy, lz, q, roots, root_class, reason)
    zm, zp = polar_roots(a, energy, lz, q)
    return KerrGeoPlunge(
        root_class,
        (a=a,),
        (E=energy, Lz=lz, Q=q),
        "Mino",
        (radial=roots, polar=(zm=zm, zp=zp)),
        (t0=0.0, radial=0.0, theta=0.0, phi0=0.0),
        nothing,
        nothing,
        (radial=(r -> _radial_potential(a, energy, lz, q, r)),
         polar_z=(z -> _polar_potential_z(a, energy, lz, q, z))),
        (radial=nothing, polar_z=nothing),
        (supported=false, reason=reason)
    )
end

"""
    kerr_geo_plunge(a, E, Lz, Q; kwargs...)

Construct a bound plunge trajectory from Kerr spin and constants of
motion. Supported root classes currently return callable trajectory fields
`t`, `r`, `theta`, `phi`, `rstar`, `u`, and `v`; status metadata records branch
support, horizon cutoff conventions, and time-origin policy.

Common keyword controls include `initPhases`, `radial_start`,
`initial_radius`, `initial_theta`, `time_origin`, and
`horizon_time_origin_offset`. Use `radial_start=:turning_point` to start at the
bound plunge turning point outside the horizon.
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
    if !(root_class in ("Complex", "Real1", "Real2"))
        return _unsupported_plunge(
            a, energy, lz, q, roots, root_class,
            "The current Julia plunge implementation supports Complex, Real1, and guarded Real2 trajectory root structures only."
        )
    end

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

    t_raw, r, theta, phi = generic_plunge_orbit(
        a, energy, lz, q;
        initPhases=phases,
        real2_horizon_offset=real2_horizon_offset,
    )
    lambda_radial_endpoint, lambda_horizon_from_turning_point, lambda_of_radius = lambda_of_r(a, energy, lz, q)
    rplus = 1 + sqrt(1 - a^2)
    horizon_lambda_from_start = lambda_horizon_from_turning_point - lambda_r0
    radial_endpoint_lambda_from_start = lambda_radial_endpoint - lambda_r0
    duration_metadata = (
        start_lambda=0.0,
        horizon_lambda=horizon_lambda_from_start,
        mino_time_to_horizon=horizon_lambda_from_start,
        radial_endpoint_lambda=radial_endpoint_lambda_from_start,
        mino_time_to_radial_endpoint=radial_endpoint_lambda_from_start,
        horizon_radius=rplus,
        coordinate_time_to_horizon_status="boyer_lindquist_t_diverges_at_future_horizon",
        retarded_time_to_horizon_status="u_diverges_linearly_at_future_horizon",
        advanced_time_to_horizon_status="v_has_finite_horizon_anchor",
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

    radial_potential(rvalue) = _radial_potential(a, energy, lz, q, rvalue)
    polar_potential_z(zvalue) = _polar_potential_z(a, energy, lz, q, zvalue)
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

    status = if root_class == "Real2"
        (
            supported=true,
            reason="ok_guarded_real2_trajectory_only",
            trajectory_support_status="guarded_real2_bound_plunge_exterior_to_horizon_cutoff",
            production_branch_support="guarded_trajectory_only_source_amplitude_waveform_blocked",
            time_phi_status="stable_realroot_analytic_branch_continuation_with_linear_anchor",
            mathematica_reference_status="public_wrapper_blocked_formula_subset_mathematica_benchmark_pass",
            source_route_status="blocked_no_real2_source_validation",
            waveform_handoff_status="metadata_only_amplitude_not_run",
            retarded_advanced_time_status="trajectory_callables_available_finite_cutoff_convention_required_near_horizon",
            real2_horizon_offset=real2_horizon_offset,
            real2_rplus=real2_rplus,
            real2_r_end=real2_r_end,
            real2_lambda_end=real2_lambda_end,
            rstar_convention=near_horizon_metadata.rstar_convention,
            near_horizon_branch=near_horizon_metadata.near_horizon_branch,
            P_plus_sign=near_horizon_metadata.P_plus_sign,
            near_horizon_series_order=near_horizon_metadata.near_horizon_series_order,
            near_horizon_series_variable=near_horizon_metadata.near_horizon_series_variable,
            u_status=near_horizon_metadata.u_status,
            v_status=near_horizon_metadata.v_status,
            lambda_to_rstar_path=near_horizon_metadata.lambda_to_rstar_path,
            series_switch_tolerance_target=near_horizon_metadata.series_switch_tolerance_target,
            series_switch_status=near_horizon_metadata.series_switch_status,
            regular_time_coefficients_status=near_horizon_metadata.regular_time_coefficients_status,
            time_origin=time_origin,
            horizon_time_origin_status=time_origin == :future_horizon_v_zero ? "future_horizon_v_zero_aligned" : "input_t0_preserved",
            horizon_time_shift=horizon_time_shift,
            horizon_v_anchor_before_shift=time_origin_anchor.vH,
            horizon_v_anchor_method=time_origin_anchor.method,
            horizon_time_origin_offset=horizon_time_origin_offset,
            duration=duration_metadata,
            horizon_lambda=horizon_lambda_from_start,
            mino_time_to_horizon=horizon_lambda_from_start,
            radial_endpoint_lambda=radial_endpoint_lambda_from_start,
            mino_time_to_radial_endpoint=radial_endpoint_lambda_from_start,
        )
    else
        (
            supported=true,
            reason="ok",
            trajectory_support_status="production_trajectory_supported",
            production_branch_support=root_class == "Complex" ? "supported_complex_current_route" : "supported_real1_orbit_trajectory_source_route_guarded",
            time_phi_status="analytic_branch_formula",
            mathematica_reference_status=root_class == "Complex" ? "available_complex_reference_record" : "available_real1_reference_record",
            source_route_status=root_class == "Complex" ? "supported_complex_current_route" : "blocked_real1_non_complex_source_route_not_validated",
            waveform_handoff_status="metadata_only_unless_source_route_validated",
            retarded_advanced_time_status="trajectory_callables_available_finite_cutoff_convention_required_near_horizon",
            real2_horizon_offset=NaN,
            real2_rplus=NaN,
            real2_r_end=NaN,
            real2_lambda_end=NaN,
            rstar_convention=near_horizon_metadata.rstar_convention,
            near_horizon_branch=near_horizon_metadata.near_horizon_branch,
            P_plus_sign=near_horizon_metadata.P_plus_sign,
            near_horizon_series_order=near_horizon_metadata.near_horizon_series_order,
            near_horizon_series_variable=near_horizon_metadata.near_horizon_series_variable,
            u_status=near_horizon_metadata.u_status,
            v_status=near_horizon_metadata.v_status,
            lambda_to_rstar_path=near_horizon_metadata.lambda_to_rstar_path,
            series_switch_tolerance_target=near_horizon_metadata.series_switch_tolerance_target,
            series_switch_status=near_horizon_metadata.series_switch_status,
            regular_time_coefficients_status=near_horizon_metadata.regular_time_coefficients_status,
            time_origin=time_origin,
            horizon_time_origin_status=time_origin == :future_horizon_v_zero ? "future_horizon_v_zero_aligned" : "input_t0_preserved",
            horizon_time_shift=horizon_time_shift,
            horizon_v_anchor_before_shift=time_origin_anchor.vH,
            horizon_v_anchor_method=time_origin_anchor.method,
            horizon_time_origin_offset=horizon_time_origin_offset,
            duration=duration_metadata,
            horizon_lambda=horizon_lambda_from_start,
            mino_time_to_horizon=horizon_lambda_from_start,
            radial_endpoint_lambda=radial_endpoint_lambda_from_start,
            mino_time_to_radial_endpoint=radial_endpoint_lambda_from_start,
        )
    end

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

"""
    kerr_geodesic(...)

Unified project-facing constructor. Use APEX-like parameters
`kerr_geodesic(a,p,e,x)` for stable-bound orbit families, or constants of
motion `kerr_geodesic(a,(E,Lz,Q))` / `kerr_geodesic(a; constants=(...))` for
bound plunge families.
"""
function kerr_geodesic(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geodesic(a, (constants.E, constants.Lz, constants.Q); kwargs...)
end

function kerr_geodesic(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    energy, lz, q = constants
    plunge = kerr_geo_plunge(a, energy, lz, q; kwargs...)
    return KerrGeodesicFamily(
        :constants,
        (a=a,),
        (E=energy, Lz=lz, Q=q),
        plunge.OrbitClass,
        nothing,
        plunge,
        (supported=plunge.Status.supported, reason=plunge.Status.reason)
    )
end

function kerr_geodesic(a::Real, p::Real, e::Real, x::Real; input::Symbol=:apex, kwargs...)
    if input == :constants
        return kerr_geodesic(a, (p, e, x); kwargs...)
    elseif input != :apex
        error("Unknown input type. Use :apex or :constants.")
    end

    constants = kerr_geo_constants_of_motion(a, p, e, x)
    energy = constants["E"]
    lz = constants["Lz"]
    q = constants["Q"]
    roots, root_class = classify_orbit(a, energy, lz, q)

    stable = nothing
    try
        if kerr_geo_orbit_type(a, p, e, x)[1] == "Bound"
            stable = kerr_geo_stable(a, p, e, x)
        end
    catch
        stable = nothing
    end

    plunge = kerr_geo_plunge(a, energy, lz, q; kwargs...)
    return KerrGeodesicFamily(
        :apex,
        (a=a, p=p, e=e, x=x),
        (E=energy, Lz=lz, Q=q),
        root_class,
        stable,
        plunge,
        (supported=(stable !== nothing || plunge.Status.supported),
        reason=plunge.Status.supported ? "ok" : plunge.Status.reason)
    )
end

function kerr_geodesic(a::Real; constants=nothing, kwargs...)
    constants === nothing && error("Provide constants=(E,Lz,Q) for one-argument kerr_geodesic.")
    return kerr_geodesic(a, constants; kwargs...)
end



end
