# The member type: one geodesic of a KerrGeodesicFamily, whatever its class. The class is the
# type parameter; the six class names below are the types users see.

"""
    KerrGeoComponent{C}

One member of a `KerrGeodesicFamily`, a single timelike geodesic in Mino time λ. `C` is its
broad class (`:stable`, `:critical`, `:plunge`, `:capture`, `:scatter`, `:trapped`);
`KerrGeoStableComponent`, `KerrGeoCriticalComponent`, … name the six types.

- `CaseId`: the case, such as `:A1`, `:K3` or `:B5`, or the tier member, such as `:A_H1` or
  `:B_X1`.
- `Tier`: `:primary`, `:horizon` (P(r₊) = 0) or `:extremal` (|a| = 1).
- `Role`: for Critical members `:on_root`, `:outer` or `:inner` (the side of the repeated root
  the member lies on); `:none` otherwise.
- `Component`: the classified radial component (`KerrGeoRadialComponent`), or `nothing` for
  horizon- and extremal-tier members and Trapped members.
- `ConstantsOfMotion`: `(a, E, Lz, Q)`.
- `Roots`: the radial and polar roots.
- `ReferenceZero`: the origins of the coordinates: the events of λ = 0 (`lambda0_event`), of
  t = φ = 0 (`t_phi_zero_event`, at `t_phi_zero_lambda` and radius `t_phi_zero_radius`) and
  of τ = 0 (`tau_zero_event`, always at λ = 0); `lambda_regular`, the λ where v = ψ = 0
  (`nothing` without a horizon-regular chart); the polar phase (`polar_phase`, or `phases`
  for Stable members).
- `Domain`: the Mino-time domain `mino` with `endpoint_closed` and `endpoint_roles`, and
  `horizon_lambda` for members that end on the future horizon.
- `Trajectory`: functions of λ: `t, r, theta, z, phi, tau`; `rstar, v, psi` where a
  horizon-regular chart exists (`v, psi` for Trapped members); `u, chi` for Trapped and
  |a| = 1 members; member-specific extras such as `lambda_of_radius` and radial increments.
- `Velocity`: the Mino-time rates `ut, ur, uz, utheta, uphi, dtau_dlambda`.
- `Potentials`, `Residuals`: R(r), Θ(z) and the geodesic-equation residuals.
- `Status`: `supported`, `spectral` (a `SpectralStatus`: the accuracy `achieved` by the
  member's Chebyshev tables and their number of `pieces`) and member-specific metadata
  (polar solution, formula family, stability; `apex`, `frequencies` and `precision` for
  Stable members).
"""
struct KerrGeoComponent{C}
    CaseId::Symbol
    Tier::Symbol
    Role::Symbol
    Component::Union{Nothing,KerrGeoRadialComponent}
    ConstantsOfMotion::NamedTuple
    Roots::NamedTuple
    ReferenceZero::NamedTuple
    Domain::NamedTuple
    Trajectory::NamedTuple
    Velocity::NamedTuple
    Potentials::NamedTuple
    Residuals::NamedTuple
    Status::NamedTuple
end

const KerrGeoStableComponent = KerrGeoComponent{:stable}
const KerrGeoCriticalComponent = KerrGeoComponent{:critical}
const KerrGeoPlungeComponent = KerrGeoComponent{:plunge}
const KerrGeoCaptureComponent = KerrGeoComponent{:capture}
const KerrGeoScatterComponent = KerrGeoComponent{:scatter}
const KerrGeoTrappedComponent = KerrGeoComponent{:trapped}

"""
    kerr_geo_member_class(m)

Broad class of member `m`, the type parameter `C` of `KerrGeoComponent{C}`: `:stable`,
`:critical`, `:plunge`, `:capture`, `:scatter` or `:trapped`.
"""
kerr_geo_member_class(::KerrGeoComponent{C}) where {C} = C

"""
    kerr_geo_sample(m, λs)

The trajectory and Mino-time four-velocity of member `m` at the Mino times `λs`: a NamedTuple
of vectors `lambda, t, r, theta, phi, tau, ut, ur, utheta, uphi`. Faster than calling the
member's functions point by point from untyped code: the loop runs behind a function barrier
on the member's concrete closures.
"""
kerr_geo_sample(m::KerrGeoComponent, λs) =
    _sample(m.Trajectory, m.Velocity, collect(Float64, λs))

function _sample(tr, u, λs::Vector{Float64})
    out = (lambda=λs, t=similar(λs), r=similar(λs), theta=similar(λs), phi=similar(λs),
        tau=similar(λs), ut=similar(λs), ur=similar(λs), utheta=similar(λs), uphi=similar(λs))
    for (i, λ) in pairs(λs)
        out.t[i] = tr.t(λ); out.r[i] = tr.r(λ); out.theta[i] = tr.theta(λ)
        out.phi[i] = tr.phi(λ); out.tau[i] = tr.tau(λ)
        out.ut[i] = u.ut(λ); out.ur[i] = u.ur(λ); out.utheta[i] = u.utheta(λ)
        out.uphi[i] = u.uphi(λ)
    end
    return out
end

# Build a member of class `class`. The tier follows the ID unless the builder knows better
# (the |a| = 1 limits of the primary cases carry primary IDs on the extremal tier).
function _member(class::Symbol, case_id::Symbol; tier=kerr_geo_tier(case_id),
        component=nothing, constants, roots=(;), reference=(;), domain,
        trajectory, velocity=(;), potentials=(;), residuals=(;), status,
        spectral=_NO_TABLES)
    role = class === :critical ? kerr_geo_critical_role(case_id) : :none
    return KerrGeoComponent{class}(case_id, tier, role, component, constants, roots,
        reference, domain, trajectory, velocity, potentials, residuals,
        merge(status, (spectral=spectral,)))
end

function Base.show(io::IO, ::MIME"text/plain", m::KerrGeoComponent{C}) where {C}
    println(io, "KerrGeo", kerr_geo_class(C).name, "Component(")
    for field in (:CaseId, :Tier, :Role, :ConstantsOfMotion, :ReferenceZero, :Domain)
        print(io, "    ", field, " = "); show(io, getfield(m, field)); println(io, ",")
    end
    print(io, "    Trajectory = ", keys(m.Trajectory), ",\n")
    print(io, "    Status = "); show(io, m.Status); println(io)
    print(io, ")")
end
