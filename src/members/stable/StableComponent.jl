# KerrGeoStable (the APEX record), KerrGeoStableComponent and the constructor `kerr_geo_stable_component`.


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

function Base.show(io::IO, kg::KerrGeoStable)
    print(io, "KerrGeoStable(constants=")
    show(io, kg.ConstantsOfMotion)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoStable)
    println(io, "KerrGeoStable")
    _show_summary_field(io, "Parameters", kg.OrbitalParameters)
    _show_summary_field(io, "Constants", kg.ConstantsOfMotion)
    _show_summary_field(io, "Orbit type", kg.OrbitalType)
    print(io, "  Trajectory = (t(lambda), r(lambda), theta(lambda), phi(lambda))")
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
    return _stable_component(a, energy, lz, q, classification, case_id; initPhases=initPhases)
end

# the Stable member of classified constants
function _stable_component(a, energy, lz, q, classification, case_id; initPhases)
    component = _stable_radial_component(classification, case_id)
    orbit, info = _class_a_orbit(a, energy, lz, q, component; initPhases=initPhases)
    stability = kerr_geo_stability_metadata(a, energy, lz, q, component)
    apex = info.apex
    residual = try
        _constants_residual((E=float(energy), Lz=float(lz), Q=float(q)),
            _apex_constants(apex.a, apex.p, apex.e, apex.x))
    catch err
        # The optional APEX reconstruction can leave its domain near E = 1.
        err isa DomainError || rethrow()
        nothing
    end
    f = info.functions
    kin = _kinematics(a, energy, lz, q, f.r, f.r, f.position, f.sign_r,
        _radial_potential_from_roots(a, energy, lz, q, classification.Status.root_structure))
    (; velocity, potentials, residuals) = _kinematic_fields(kin)
    return _member(:stable, component.CaseId, kerr_geo_tier(component.CaseId), component,
        (a=float(a), E=float(energy), Lz=float(lz), Q=float(q)),
        (radial=info.roots, polar=info.polar),
        # λ = 0 is the event of the initial phases (periapsis and the northern polar turning
        # point for zero phases); t(0) = qt0, φ(0) = qφ0, τ(0) = 0
        (lambda0_event=:initial_phases, t_phi_zero_event=:initial_phases,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=f.r(0.0), tau_zero_event=:initial_phases,
            lambda_regular=nothing, phases=orbit.InitialPhases),
        (mino=(-Inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:infinite_past_worldline, :infinite_future_worldline)),
        (t=f.t, r=f.r, theta=f.theta, z=f.z, phi=f.phi, tau=f.tau),
        velocity, potentials, residuals,
        (supported=true, formula_family=component.FormulaFamily,
            formula_kind=component.CaseId === :A1 ? :jacobi_libration : :constant_radius,
            shape=stability.shape, stability=stability.stability,
            stability_check_passed=stability.stability_check_passed, limit=stability.limit,
            radial_derivatives=stability.radial_derivatives, apex=apex,
            frequencies=orbit.Frequencies, orbit=orbit, constants_residual=residual,
            polar=info.polar, precision=_stable_precision(orbit.Frequencies)),
        info.spectral)
end

# t grows by ϒt·2π/ϒr per radial period, so after one period t is known to one ulp of that
# at best: ~1e3 M for |E − 1| ~ 1e-13 (apoapsis ~1e13 M), whatever the method
function _stable_precision(frequencies)
    period_t = frequencies.ϒt * (2 * oftype(frequencies.ϒt, π)) / frequencies.ϒr
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
