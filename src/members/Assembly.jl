# The assembly shared by the members built on the radial coordinate engine. A member supplies
# its radial track and polar solution; z, θ come from the polar solution, t, φ, τ (and r*, v, ψ
# where the member has a horizon-regular chart) from the engine, the Mino-time velocity and the
# geodesic-equation residuals from `_kinematics`.
#
# A track is a NamedTuple with
#   r(λ), check(λ), check_bl(λ)   radius and the Mino-time domains of (r, z, τ, v, ψ) and (t, φ)
#   coords                        `_engine_coordinates(...)`, or (t, phi, tau) closures
#   sign_r(λ)                     sign of dr/dλ
#   domain, reference             Domain and ReferenceZero (reference.lambda_regular is the λ
#                                 where v = ψ = 0, or `nothing` without a regular chart)
#   trajectory                    member-specific extras: λ(r), radial increments, …

function _constant_radial_metadata(radius)
    return (domain=(mino=(-Inf, Inf), endpoint_closed=(false, false),
            endpoint_roles=(:infinite_past_worldline, :infinite_future_worldline)),
        reference=(lambda0_event=:polar_phase_reference, t_phi_zero_event=:polar_phase_reference,
            t_phi_zero_lambda=0.0, t_phi_zero_radius=radius, tau_zero_event=:polar_phase_reference,
            lambda_regular=nothing))
end

# Preserve NamedTuple merge order without compiling a merge of nested closure trees.
Base.@nospecializeinfer @noinline function _merge_member_fields(@nospecialize(parts::Tuple))
    fields = Symbol[]
    items = Any[]
    for part in parts, name in fieldnames(typeof(part))
        index = findfirst(==(name), fields)
        value = getfield(part, name)
        if index === nothing
            push!(fields, name)
            push!(items, value)
        else
            items[index] = value
        end
    end
    return NamedTuple{Tuple(fields)}(Tuple(items))
end

# The coordinate functions of a track: an `EngineCoordinates` gives one closure per coordinate,
# each capturing only that object; a track that overrides coordinates passes a NamedTuple of
# functions (at least t, phi, tau).
_coordinate_functions(c::EngineCoordinates) = (t=λ -> _coords_t(c, λ), phi=λ -> _coords_phi(c, λ),
    tau=λ -> _coords_tau(c, λ), v=λ -> _coords_v(c, λ), psi=λ -> _coords_psi(c, λ),
    radial=λ -> _coords_radial(c, λ), spectral=() -> _coords_spectral(c))
_coordinate_functions(c::NamedTuple) = c
_coordinate_rstar(c::NamedTuple,a,rbl,lambda) = kerr_rstar(a,rbl(lambda))
function _coordinate_rstar(c::EngineCoordinates,a,rbl,lambda)
    radius=rbl(lambda)
    state=_radial_state(c.r_of,lambda)
    return state === nothing ? kerr_rstar(a,radius) : _relative_rstar(state)
end

"""
    _engine_member(class, id, a, E, Lz, Q, polar, track, component, structure, potential, tier,
                   roots, status, chart)

Assemble a member of class `class` from its radial `track` and its `polar` solution;
`structure` is the root structure (for the product-form potential when `potential` is
`nothing`), `chart` the unshifted (r_*, φ_H) of members without a regular chart, or `nothing`.
Compiled once: nothing numerical happens here, and the closures it creates keep the concrete
types of what they capture.
"""
Base.@nospecializeinfer @noinline function _engine_member(class, id, a, energy, lz, q,
        @nospecialize(polar), @nospecialize(track), @nospecialize(component),
        @nospecialize(structure), @nospecialize(potential), tier, @nospecialize(roots),
        @nospecialize(status), @nospecialize(chart))
    r, check, check_bl = getfield(track, :r), getfield(track, :check), getfield(track, :check_bl)
    cf = _coordinate_functions(getfield(track, :coords))
    tf, phif, tauf = getfield(cf, :t), getfield(cf, :phi), getfield(cf, :tau)
    rbl(λ) = r(check_bl(λ))
    position = _polar_position(polar)
    z(λ) = position(check(λ))[1]
    reference = getfield(track, :reference)
    regular = reference.lambda_regular !== nothing
    t(λ) = tf(check_bl(λ))
    phi(λ) = phif(check_bl(λ))
    # horizon charts: from the engine when the member reaches a horizon, otherwise the unshifted
    # v = t + r_*, ψ = φ + φ_H (and u = t − r_*, χ = φ − φ_H) when a `chart` (r_*, φ_H) is given
    charts = if regular
        vf, psif = getfield(cf, :v), getfield(cf, :psi)
        coords=getfield(track,:coords)
        (rstar=λ -> _coordinate_rstar(coords,a,rbl,λ), v=λ -> vf(check(λ)), psi=λ -> psif(check(λ)))
    elseif chart === nothing
        (;)
    else
        rstarf, phihf = chart.rstar, chart.phi_h
        (rstar=λ -> rstarf(rbl(λ)), v=λ -> t(λ) + rstarf(rbl(λ)),
            psi=λ -> phi(λ) + phihf(rbl(λ)), u=λ -> t(λ) - rstarf(rbl(λ)),
            chi=λ -> phi(λ) - phihf(rbl(λ)))
    end
    radial = if haskey(cf, :radial)
        radialf = getfield(cf, :radial)
        (radial_t=λ -> radialf(check_bl(λ))[1], radial_phi=λ -> radialf(check_bl(λ))[2],
            radial_tau=λ -> radialf(check(λ))[3])
    else
        (;)
    end
    trajectory = _merge_member_fields(((
        t=t,
        r=r,
        theta=λ -> (p=position(check(λ)); atan(sqrt(p[3]), p[1])),
        z=z,
        phi=phi,
        tau=λ -> tauf(check(λ)),
        ), charts, radial, getfield(track, :trajectory),
    ))
    kin = _kinematics(a, energy, lz, q, r, rbl, λ -> position(check(λ)), getfield(track, :sign_r),
        potential === nothing ? _radial_potential_from_roots(a, energy, lz, q, structure) : potential;
        radial_track=haskey(track,:radial_track) ? getfield(track,:radial_track) : nothing)
    (; velocity, potentials, residuals) = _kinematic_fields(kin)
    spectralf = haskey(cf, :spectral) ? getfield(cf, :spectral) : () -> (achieved=0.0, pieces=0)
    return _member(class, id, tier, component, (a=float(a), E=float(energy), Lz=float(lz),
            Q=float(q)),
        merge(component === nothing ? (;) :
            (radial=Tuple(item.radius for item in component.Metadata.roots),),
            roots, (polar=polar.metadata,)), reference, getfield(track, :domain), trajectory,
        velocity, potentials, residuals, merge((supported=true,), status, (polar=polar.metadata,)),
        SpectralStatus(() -> (spectralf(), _polar_spectral(polar))))
end

"""The classified component of class `class` of these constants (the case `requested`, if given)."""
function _class_component(classification, class, requested=nothing)
    candidates = [c for c in classification.Components
        if c.BroadClass === class && c.CaseId !== nothing &&
            (requested === nothing || c.CaseId === requested)]
    name = kerr_geo_class(class).name
    isempty(candidates) && error("These constants have no $(name) component" *
        (requested === nothing ? "" : " $(requested)") * "; their cases are " *
        "$(classification.CaseIds).")
    length(candidates) == 1 || error("The $(name) component is ambiguous: " *
        "$([c.CaseId for c in candidates]); select one with case_id.")
    return only(candidates)
end

# the Velocity, Potentials and Residuals fields of a member from its `_kinematics`
_kinematic_fields(kin) = (
    velocity=(ut=kin.ut, ur=kin.ur, uz=kin.uz, utheta=kin.utheta, uphi=kin.uphi,
        dtau_dlambda=kin.dtau_dlambda),
    potentials=(radial=kin.R, polar_z=kin.Θ),
    residuals=(radial=kin.radial_residual, polar_z=kin.polar_residual,
        normalization=kin.normalization_residual),
)
