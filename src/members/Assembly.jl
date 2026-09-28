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

"""Assemble a member of class `class` from its radial `track` and its `polar` solution."""
function _engine_member(class, id, a, energy, lz, q, polar, track; component=nothing,
        tier=kerr_geo_tier(id), roots=(;), status=(;))
    (; r, check, check_bl, coords) = track
    rbl(λ) = r(check_bl(λ))
    position = _polar_position(polar)
    z(λ) = position(check(λ))[1]
    regular = track.reference.lambda_regular !== nothing
    trajectory = (
        t=λ -> coords.t(check_bl(λ)),
        r=r,
        theta=λ -> acos(clamp(z(λ), -1.0, 1.0)),
        z=z,
        phi=λ -> coords.phi(check_bl(λ)),
        tau=λ -> coords.tau(check(λ)),
        (regular ? (rstar=λ -> kerr_rstar(a, rbl(λ)), v=λ -> coords.v(check(λ)),
            psi=λ -> coords.psi(check(λ))) : (;))...,
        (haskey(coords, :radial) ? (radial_t=λ -> coords.radial(check_bl(λ))[1],
            radial_phi=λ -> coords.radial(check_bl(λ))[2],
            radial_tau=λ -> coords.radial(check(λ))[3]) : (;))...,
        track.trajectory...,
    )
    kin = _kinematics(a, energy, lz, q; r=r, rbl=rbl, z=z, uz=λ -> position(check(λ))[2],
        sin2=λ -> position(check(λ))[3], sign_r=track.sign_r, R=rv -> kerr_radial_potential(a, energy, lz, q, rv))
    (; velocity, potentials, residuals) = _kinematic_fields(kin)
    return _member(class, id; tier=tier, component=component,
        constants=(a=float(a), E=float(energy), Lz=float(lz), Q=float(q)),
        roots=merge(component === nothing ? (;) :
            (radial=Tuple(item.radius for item in component.Metadata.roots),),
            roots, (polar=polar.metadata,)), reference=track.reference,
        domain=track.domain, trajectory, velocity, potentials, residuals,
        status=merge((supported=true,), status, (polar=polar.metadata,)),
        spectral=SpectralStatus(() -> (haskey(coords, :spectral) ? coords.spectral() :
            (achieved=0.0, pieces=0), _polar_spectral(polar))))
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

# the Velocity, Potentials and Residuals fields of a member from its `_kinematics`. Members
# pass them to `_member` as three keywords, not by splatting: a splatted keyword list is
# merged into the other keywords, and `merge` of NamedTuples this large (closures capturing the
# coordinate engine) takes minutes to compile.
_kinematic_fields(kin) = (
    velocity=(ut=kin.ut, ur=kin.ur, uz=kin.uz, utheta=kin.utheta, uphi=kin.uphi,
        dtau_dlambda=kin.dtau_dlambda),
    potentials=(radial=kin.R, polar_z=kin.Θ),
    residuals=(radial=kin.radial_residual, polar_z=kin.polar_residual,
        normalization=kin.normalization_residual),
)
