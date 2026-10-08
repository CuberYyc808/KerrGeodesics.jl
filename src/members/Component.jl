# The member type: one geodesic of a KerrGeodesicFamily, whatever its class. The class is the
# type parameter; the six class names below are the types users see.

"""
    KerrGeoComponent{C}

One member of a `KerrGeodesicFamily`, a single timelike geodesic in Mino time ``\\lambda``. `C` is its
broad class (`:stable`, `:critical`, `:plunge`, `:capture`, `:scatter`, `:trapped`);
`KerrGeoStableComponent`, `KerrGeoCriticalComponent`, … name the six types.

- `CaseId`: the case, such as `:A1`, `:K3` or `:B5`, or the tier member, such as `:A_H1` or
  `:B_X1`.
- `Tier`: `:primary`, `:horizon` (``P(r_+) = 0``) or `:extremal` (``|a| = 1``).
- `Role`: for Critical members `:on_root`, `:outer` or `:inner` (the side of the repeated root
  the member lies on); `:none` otherwise.
- `Component`: the classified radial component (`KerrGeoRadialComponent`), or `nothing` for
  horizon- and extremal-tier members and Trapped members.
- `ConstantsOfMotion`: `(a, E, Lz, Q)`.
- `Roots`: the radial and polar roots.
- `ReferenceZero`: the origins of the coordinates: the events of ``\\lambda = 0`` (`lambda0_event`), of
  ``t = \\phi = 0`` (`t_phi_zero_event`, at `t_phi_zero_lambda` and radius `t_phi_zero_radius`) and
  of ``\\tau = 0`` (`tau_zero_event`, always at ``\\lambda = 0``); `lambda_regular`, the ``\\lambda`` where ``v = \\psi = 0``
  (`nothing` without a horizon-regular chart); the polar phase (`polar_phase`, or `phases`
  for Stable members).
- `Domain`: the Mino-time domain `mino` with `endpoint_closed` and `endpoint_roles`, and
  `horizon_lambda` for members that end on the future horizon.
- `Trajectory`: functions of ``\\lambda``: `t, r, theta, z, phi, tau` (``z = \\cos \\theta``); `rstar, v, psi`
  (``r_*``, ``v = t + r_*``, ``\\psi = \\phi + \\phi_H``) where a horizon-regular chart exists (`v, psi` for Trapped members); `u, chi` for Trapped and
  ``|a| = 1`` members; member-specific extras such as `lambda_of_radius` and radial increments.
- `Velocity`: the Mino-time rates, functions of ``\\lambda``: `ut` ``= dt/d\\lambda``, `ur` ``= dr/d\\lambda = \\pm\\sqrt{R(r)}``,
  `uz` ``= dz/d\\lambda = \\pm\\sqrt{\\Theta(z)}``, `utheta` ``= d\\theta/d\\lambda = -(dz/d\\lambda)/\\sin \\theta``, `uphi` ``= d\\phi/d\\lambda`` and
  `dtau_dlambda` ``= d\\tau/d\\lambda = \\Sigma = r^2 + a^2 z^2``; the four-velocity is ``u^\\mu = (dx^\\mu/d\\lambda)/\\Sigma``.
- `Potentials`: ``R(r)`` (`radial`) and ``\\Theta(z)`` (`polar_z`) as functions of their variable.
- `Residuals`: functions of ``\\lambda``: ``(dr/d\\lambda)^2 - R``, ``(dz/d\\lambda)^2 - \\Theta`` and ``g_{\\mu\\nu}u^\\mu u^\\nu + 1``.
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
    # the record fields are abstract: one constructor for every argument type (the default one
    # would be compiled for each member's closure types)
    Base.@nospecializeinfer function KerrGeoComponent{C}(case_id, tier, role, @nospecialize(component),
            @nospecialize(constants), @nospecialize(roots), @nospecialize(reference),
            @nospecialize(domain), @nospecialize(trajectory), @nospecialize(velocity),
            @nospecialize(potentials), @nospecialize(residuals), @nospecialize(status)) where {C}
        return new{C}(case_id, tier, role, component, constants, roots, reference, domain,
            trajectory, velocity, potentials, residuals, status)
    end
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
of vectors `lambda, t, r, theta, phi, tau, ut, ur, utheta, uphi`, where `ut` ``= dt/d\\lambda``,
`ur` ``= dr/d\\lambda``, `utheta` ``= d\\theta/d\\lambda`` and `uphi` ``= d\\phi/d\\lambda`` (divide by
``\\Sigma = r^2 + a^2\\cos^2\\theta`` for ``dx^\\mu/d\\tau``). Faster than calling the
member's functions point by point from untyped code: the loop runs behind a function barrier
on the member's concrete closures.
"""
function kerr_geo_sample(m::KerrGeoComponent, λs)
    tr, u = m.Trajectory, m.Velocity
    return _sample(collect(typeof(m.ConstantsOfMotion.E), λs), tr.t, tr.r, tr.theta, tr.phi,
        tr.tau, u.ut, u.ur, u.utheta, u.uphi)
end

# (the barrier is specialized on the nine functions it calls, not on the whole records)
function _sample(λs::Vector{<:Real}, t, r, theta, phi, tau, ut, ur, utheta, uphi)
    out = (lambda=λs, t=similar(λs), r=similar(λs), theta=similar(λs), phi=similar(λs),
        tau=similar(λs), ut=similar(λs), ur=similar(λs), utheta=similar(λs), uphi=similar(λs))
    for (i, λ) in pairs(λs)
        out.t[i] = t(λ); out.r[i] = r(λ); out.theta[i] = theta(λ)
        out.phi[i] = phi(λ); out.tau[i] = tau(λ)
        out.ut[i] = ut(λ); out.ur[i] = ur(λ); out.utheta[i] = utheta(λ)
        out.uphi[i] = uphi(λ)
    end
    return out
end

# Build a member of class `class` (positional, compiled once: the record fields of
# KerrGeoComponent are abstract, so nothing is gained by specializing on the closures). `tier`
# follows the ID unless the builder knows better (the |a| = 1 limits of the primary cases carry
# primary IDs on the extremal tier); `spectral` is the member's `SpectralStatus`.
Base.@nospecializeinfer @noinline function _member(class::Symbol, case_id::Symbol, tier::Symbol,
        @nospecialize(component), @nospecialize(constants), @nospecialize(roots),
        @nospecialize(reference), @nospecialize(domain), @nospecialize(trajectory),
        @nospecialize(velocity), @nospecialize(potentials), @nospecialize(residuals),
        @nospecialize(status), @nospecialize(spectral))
    role = class === :critical ? kerr_geo_critical_role(case_id) : :none
    T = _float_type(values(constants)...)
    p = precision(T(constants.E))
    return KerrGeoComponent{class}(case_id, tier, role, component, _retype(T, constants),
        _retype(T, roots), _retype(T, reference), _retype(T, domain),
        _precision_wrap(T, p, trajectory), _precision_wrap(T, p, velocity),
        _precision_wrap(T, p, potentials), _precision_wrap(T, p, residuals),
        _retype(T, _merge_member_fields((status, (spectral=spectral,)))))
end

# the numbers of a member's records in its floating-point type T (literal 0.0, ±Inf and NaN
# of the builders included); integers, symbols and functions are kept
_retype(::Type{T}, x::AbstractFloat) where {T} = T(x)
_retype(::Type{T}, x::Complex{<:AbstractFloat}) where {T} = Complex{T}(x)
_retype(::Type{T}, x::Union{Tuple,NamedTuple}) where {T} = map(v -> _retype(T, v), x)
_retype(::Type{T}, x::AbstractVector{<:AbstractFloat}) where {T} = T.(x)
_retype(::Type, x) = x

function Base.show(io::IO, m::KerrGeoComponent{C}) where {C}
    print(io, "KerrGeo", kerr_geo_class(C).name, "Component(", m.CaseId,
        ", tier=", m.Tier, ", constants=")
    show(io, m.ConstantsOfMotion)
    print(io, ")")
end

function _show_summary_field(io::IO, label, value)
    print(io, "  ", rpad(label, 10), " = ")
    show(IOContext(io, :compact => true, :limit => true), value)
    println(io)
end

function _show_summary_status(io::IO, status)
    supported = get(status, :supported, nothing)
    print(io, "  Status     = ", supported === nothing ? "not recorded" :
        supported ? "supported" : "unsupported")
    errors = get(status, :member_errors, ())
    isempty(errors) || print(io, "; ", length(errors), " member error(s)")
    reason = get(status, :reason, nothing)
    supported === false && reason !== nothing && print(io, "; ", reason)
end

function Base.show(io::IO, ::MIME"text/plain", m::KerrGeoComponent{C}) where {C}
    println(io, "KerrGeo", kerr_geo_class(C).name, "Component (", m.CaseId, ")")
    _show_summary_field(io, "Constants", m.ConstantsOfMotion)
    _show_summary_field(io, "Tier", m.Tier)
    C === :critical && _show_summary_field(io, "Role", m.Role)
    _show_summary_field(io, "Mino time", get(m.Domain, :mino, nothing))
    println(io, "  Trajectory = (t(lambda), r(lambda), theta(lambda), phi(lambda))")
    _show_summary_status(io, m.Status)
end
