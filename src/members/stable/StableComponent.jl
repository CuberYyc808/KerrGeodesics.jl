# KerrGeoStable (APEX record), KerrGeoStableComponent and the constructors `kerr_geo_stable`,
# `kerr_geo_stable_component`.

"""
    KerrGeoStable

Record of a Stable (A1/A2) Kerr geodesic in APEX form: parameters `(a, p, e, x)`, constants,
Mino frequencies, and callable t, r, θ, φ, four-velocity and cross functions of Mino time.
`kerr_geo_stable` returns it; `kerr_geo_stable_component` stores it as `Status.orbit`.
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

The stable orbit with APEX parameters `(a, p, e, x)` as a `KerrGeoStable` record. The
initial phases `(qt0, qr0, qθ0, qφ0)` shift t, the radial phase, the polar phase and φ at
λ = 0; with zero phases, λ = 0 is at periapsis and at the northern polar turning point.
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

"""
    kerr_geo_stable_component(a, E, Lz, Q; case_id=nothing, initPhases=(0,0,0,0))

The Stable member (A1 eccentric, A2 circular/spherical) of the constants `(a, E, Lz, Q)`, a
`KerrGeoStableComponent` built from the turning points of the classified radial component.
With zero phases, λ = 0 is at periapsis and at the northern polar turning point;
`initPhases = (qt0, qr0, qθ0, qφ0)` shifts the phases as in `kerr_geo_orbit`, and τ(0) = 0.
`Status` carries the APEX parameters (p, e, x) derived from the turning points
(`Status.apex`), the Mino frequencies (`Status.frequencies`), the `KerrGeoStable` record
(`Status.orbit`), `Status.precision` (the BL time of one radial period and its ulp) and
`Status.constants_residual`, the difference between the input constants and those recomputed
from (p, e, x).
"""
function kerr_geo_stable_component(a::Real, energy::Real, lz::Real, q::Real;
        case_id=nothing,
        initPhases=(0.0, 0.0, 0.0, 0.0))
    classification = kerr_geo_classify(a, energy, lz, q)
    component = _stable_radial_component(classification, case_id)
    orbit, info = _class_a_orbit(a, energy, lz, q, component; initPhases=initPhases)
    stability = kerr_geo_stability_metadata(a, energy, lz, q, component)
    apex = info.apex
    residual = try
        _constants_residual((E=float(energy), Lz=float(lz), Q=float(q)),
            kerr_geo_constants_of_motion(apex.a, apex.p, apex.e, apex.x))
    catch
        nothing
    end
    f = info.functions
    kin = _kinematics(a, energy, lz, q; r=f.r, z=f.z, uz=f.uz, sin2=f.sin2, sign_r=f.sign_r,
        R=rv -> kerr_radial_potential(a, energy, lz, q, rv))
    (; velocity, potentials, residuals) = _kinematic_fields(kin)
    return _member(:stable, component.CaseId; component=component,
        constants=(a=float(a), E=float(energy), Lz=float(lz), Q=float(q)),
        roots=(radial=info.roots, polar=info.polar),
        # λ = 0 is the event of the initial phases (periapsis and the northern polar turning
        # point for zero phases); t(0) = qt0, φ(0) = qφ0, τ(0) = 0
        reference=(lambda0_event=:initial_phases, t_phi_zero_event=:initial_phases,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=f.r(0.0), tau_zero_event=:initial_phases,
            lambda_regular=nothing, phases=orbit.InitialPhases),
        domain=(mino=(-Inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:infinite_past_worldline, :infinite_future_worldline)),
        trajectory=(t=f.t, r=f.r, theta=f.theta, z=f.z, phi=f.phi, tau=f.tau),
        velocity, potentials, residuals,
        status=(supported=true, formula_family=component.FormulaFamily,
            formula_kind=component.CaseId === :A1 ? :jacobi_libration : :constant_radius,
            shape=stability.shape, stability=stability.stability,
            stability_check_passed=stability.stability_check_passed, limit=stability.limit,
            radial_derivatives=stability.radial_derivatives, apex=apex,
            frequencies=orbit.Frequencies, orbit=orbit, constants_residual=residual,
            polar=info.polar, precision=_stable_precision(orbit.Frequencies)),
        spectral=info.spectral)
end

# t grows by ϒt·2π/ϒr per radial period, so after one period t is known to one ulp of that
# at best: ~1e3 M for |E − 1| ~ 1e-13 (apoapsis ~1e13 M), whatever the method
function _stable_precision(frequencies)
    period_t = frequencies.ϒt * 2pi / frequencies.ϒr
    return (t_radial_period=period_t, t_ulp_per_period=eps(period_t))
end

function kerr_geo_stable_component(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_stable_component(
        a, constants.E, constants.Lz, constants.Q; kwargs...)
end

function kerr_geo_stable_component(
        a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_stable_component(a, constants...; kwargs...)
end
