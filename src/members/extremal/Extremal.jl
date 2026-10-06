# Exact-extremal (|a| = 1) members: classification, radial models and assembly, for the X-tier
# members A-X1, A-X2, B-X1, B-X2, C-X1..C-X4, D-X1, D-X2 and the primary cases at |a| = 1.
# Primary members use the shared coordinate engine. The P_H = 0 members retain
# their elementary horizon-root primitives; the polar parts use the polar engine.
#
# Zero conventions (`ReferenceZero`) follow the other members except where the |a| = 1
# construction anchors differently: every member that ends on the future horizon (B, C, K2,
# K5, K8, K11 at |a| = 1) has λ = 0 and v = ψ = 0 on the future horizon, where t and φ
# diverge, so t and φ are fixed by that chart (`t_phi_zero_event =
# :regular_chart_at_future_horizon`) rather than at a finite reference radius; A1 at |a| = 1
# has λ = 0 and t = φ = 0 at the inner turning point instead of APEX initial phases; members
# that do not reach a horizon use v = t + r_*, ψ = φ + φ_H with no additive shift
# (`lambda_regular = nothing`).

"""
    KerrGeoExtremalFamily

All members admitted by one set of constants at exact `a = ±1`, returned by
`kerr_geo_extremal_family`. Fields: `MetricLimit` (`:extremal_plus` or `:extremal_minus`),
`ConstantsOfMotion` `(a, E, Lz, Q)`, `Classification` (the admitted `case_ids`, the
`excluded` outcomes and `P_H`), `Members` (a tuple of components) and `Status`.
"""
struct KerrGeoExtremalFamily
    MetricLimit::Symbol
    ConstantsOfMotion::NamedTuple
    Classification::NamedTuple
    Members::Tuple
    Status::NamedTuple
end

function Base.show(io::IO, family::KerrGeoExtremalFamily)
    print(io, "KerrGeoExtremalFamily(", family.MetricLimit, ", cases=")
    show(io, family.Classification.case_ids)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", family::KerrGeoExtremalFamily)
    println(io, "KerrGeoExtremalFamily (", family.MetricLimit, ")")
    _show_summary_field(io, "Constants", family.ConstantsOfMotion)
    _show_summary_field(io, "Cases", family.Classification.case_ids)
    _show_summary_field(io, "Members", length(family.Members))
    _show_summary_status(io, family.Status)
end


# The horizon-root radial motion in z = 1/(r − 1): dλ = −dz/√h with the quadratic
# h(z) = h2 z² + h1 z + h0, D = √(h1² − 4 h2 h0), and its primitives k_j = ∫ z^j dz/√h. With
# w = √h and the position a = (2 h2 z + h1)/D on the parabola (|a| = 1 at a root of h, a = 0 at
# its vertex), one primitive serves every h2:
#   a > 1/√2 (the branch that meets a root at a = 1, and all of h2 > 0): written without 1/h2,
#   so it holds uniformly as h2 → 0 (the linear limit) — ρ = D + h1, X = 4 h2 w²/D²,
#     k0 = 2 (w/D) S(X), k1 = −4 h0 w/(D ρ) − 4 h1 w³ S₁(X)/D³,
#     k2 = w [−2 h0 z/(D ρ) − z²/(2D) + 6 h0²/(D ρ²)] − 4 h0 w³/(3 D³) + 4 w⁵ (3 D² + 8 h2 h0) S₂(X)/D⁵,
#     S(X) = asinh(√X)/√X (asin(√−X)/√−X for X < 0), S₁ = (S − 1)/X, S₂ = (S₁ + 1/6)/X;
#   a ≤ 1/√2 (h2 < 0 only: the arc of the concave parabola through its vertex, where the
#   vertex is far from the linear limit): the same function k0 = acos(a)/√(−h2) and the
#   classical k1 = w/h2 − h1 k0/(2 h2), k2 from ∫ d(z w).
# Both pieces are the one primitive (acos a = asin(√(1 − a²)) for a > 0), so differences
# across the switch carry no constant. The same functions serve the reversed quadratic
# h0 u² + h1 u + h2 in u = 1/z. At a turning point (a root of h) w is exactly zero and |a|
# exactly one, so the primitives take their branch values without the √eps error of a
# rounded root.
function _horizon_root_S(X)
    T = float(typeof(X))
    if abs(X) < 0.25
        # asinh(√X)/√X = Σ c_k X^k, c_{k+1}/c_k = −(2k + 1)²/((2k + 2)(2k + 3)); S₁ and S₂ from the
        # same terms. Each term is at most 1/4 of the previous: 33 terms (tail below 4^(−30)) in
        # Float64, ⌈p/2⌉ + 6 for p bits
        c = one(T); S = zero(T); S1 = zero(T); S2 = zero(T); Xk = one(T)
        for k in 0:cld(precision(T), 2) + 5
            S += c * Xk
            k >= 1 && (S1 += c * Xk / X)
            k >= 2 && (S2 += c * Xk / X^2)
            c *= T(-(2k + 1)^2) / ((2k + 2) * (2k + 3))
            Xk *= X
        end
        return iszero(X) ? (one(T), -inv(T(6)), T(3) / 40) : (S, S1, S2)
    end
    S = X > 0 ? asinh(sqrt(X)) / sqrt(X) : asin(sqrt(-X)) / sqrt(-X)
    S1 = (S - 1) / X
    return S, S1, (S1 + inv(T(6))) / X
end

# D and ρ = D + h1 of the quadratic, ρ without cancellation for h1 < 0
function _horizon_root_scales(h2, h1, h0)
    D = sqrt(h1^2 - 4h2 * h0)
    return D, (h1 >= 0 ? D + h1 : -4h2 * h0 / (D - h1))
end

function _horizon_root_primitives(h2, h1, h0, z; root=false, sqrt_h=nothing)
    h = (h2 * z + h1) * z + h0
    # a negative h beyond the rounding of its own terms lies outside the allowed interval
    T = typeof(h)
    root || h >= -32 * eps(T) * (abs(h2 * z^2) + abs(h1 * z) + abs(h0)) || throw(DomainError(h,
        "The horizon-root quadratic lies outside its allowed radial interval."))
    w = root ? zero(T) : sqrt_h === nothing ? sqrt(max(h, zero(T))) : sqrt_h
    D, ρ = _horizon_root_scales(h2, h1, h0)
    position = root ? sign(2h2 * z + h1) : (2h2 * z + h1) / D
    if position > 1 / sqrt(T(2))
        S, S1, S2 = _horizon_root_S(4h2 * w^2 / D^2)
        k0 = 2 * (w / D) * S
        k1 = -4h0 * w / (D * ρ) - 4h1 * w^3 * S1 / D^3
        k2 = w * (-2h0 * z / (D * ρ) - z^2 / (2D) + 6h0^2 / (D * ρ^2)) - 4h0 * w^3 / (3D^3) +
            4 * w^5 * (3D^2 + 8h2 * h0) * S2 / D^5
        return k0, k1, k2
    end
    h2 < 0 || throw(DomainError(z, "The horizon-root primitives of h2 > 0 lie on the branch 2 h2 z + h1 ≥ D."))
    # acos(a): near a = −1 through 1 + a = −4 h2 h/(D (D − 2 h2 z − h1)), which keeps its digits
    angle = position < -1 / sqrt(T(2)) ?
        pi - 2 * asin(sqrt(-2h2 * w^2 / (D * (D - 2h2 * z - h1)))) : acos(position)
    k0 = angle / sqrt(-h2)
    k1 = w / h2 - h1 * k0 / (2h2)
    k2 = ((2h2 * z + h1) * w / (4h2) - D^2 * k0 / (8h2) - h1 * k1 - h0 * k0) / h2
    return k0, k1, k2
end

# (cosh √y − 1)/y and sinh(√y)/√y for y of either sign (cos, sin for y < 0), the inversion of
# k0: 2 h2 z + h1 = D cosh(√h2 k0)
_horizon_root_G(y) = y > 0 ? 2 * sinh(sqrt(y) / 2)^2 / y : y < 0 ? 2 * sin(sqrt(-y) / 2)^2 / (-y) :
    one(float(y)) / 2
_horizon_root_Sh(y) = y > 0 ? sinh(sqrt(y)) / sqrt(y) : y < 0 ? sin(sqrt(-y)) / sqrt(-y) :
    one(float(y))

# The horizon-root radial model in z = 1/(r − 1) (header of `_horizon_root_primitives`). Near the
# horizon root r = 1 + 1/z loses all of z's digits (and saturates at r = 1), so the Mino-time
# path works in z: `_hr_basis_from(m, left, delta)` evaluates the primitives at the z reached
# after delta. One struct holds the quadratic and the branch's turning points; the functions
# below evaluate it (closures over each other would nest the model in every member type).
struct _HorizonRootModel{T,U}
    h2::T; h1::T; h0::T
    D::T; ρ::T
    k0inf::T                        # k0 at z = 0 (r = ∞) when the orbit reaches infinity
    turns::U                        # the branch's turning points (roots of h)
end

function _horizon_root_model(energy,q; turns=())
    T=_float_type(energy,q)
    h0=_e2m1(energy); h1=4energy^2-2; h2=3energy^2-1-q
    D,ρ=_horizon_root_scales(h2,h1,h0)
    k0inf=h0>=0 ? _horizon_root_primitives(h2,h1,h0,zero(T))[1] : T(NaN)
    return _HorizonRootModel{T,typeof(turns)}(h2,h1,h0,D,ρ,k0inf,turns)
end

# (1/z is a root of the reversed quadratic exactly when z is a root of h)
function _hr_basis_z(m::_HorizonRootModel,z; root=false, mino=nothing, sqrt_h=nothing)
    h2,h1,h0=m.h2,m.h1,m.h0
    k0,k1,k2=_horizon_root_primitives(h2,h1,h0,z;root=root,sqrt_h=sqrt_h)
    km0,km1,_=_horizon_root_primitives(h0,h1,h2,inv(z);root=root,
        sqrt_h=sqrt_h === nothing ? nothing : sqrt_h / abs(z))
    # On a Mino-time trajectory k0 is already known; reconstructing it from
    # the rounded inverse coordinate loses digits in a thin island.
    mino===nothing || (k0=-mino)
    return (I0=-k0,I1=-(k0-km0),I2=-(k0-2km0-km1),J1=-k1,J2=-k2)
end
# the branch's turning points (the same radii its spec was built with) are roots of h
_hr_basis(m::_HorizonRootModel,r)=_hr_basis_z(m,inv(r-1); root=r in m.turns)
# z(k0) inverts k0(z): 2 h2 z + h1 = D cosh(√h2 k0), written as z = D k0² G(h2 k0²)/2 − 2 h0/ρ;
# when the orbit reaches infinity (h0 >= 0) as the difference from z(k0∞) = 0, so z -> 0
# (r -> infinity) keeps its digits
function _hr_z_of_target(m::_HorizonRootModel,target)
    h2,h0,D,ρ,k0inf=m.h2,m.h0,m.D,m.ρ,m.k0inf
    k0=-target
    h0>=0 && return D*(k0+k0inf)*(k0-k0inf)*_horizon_root_Sh(h2*(k0+k0inf)^2/4)*
        _horizon_root_Sh(h2*(k0-k0inf)^2/4)/4
    return D*k0^2*_horizon_root_G(h2*k0^2)/2-2h0/ρ
end
_hr_inverse_from(m::_HorizonRootModel,left,delta)=1+inv(_hr_z_of_target(m,_hr_basis(m,left).I0+delta))
function _hr_basis_from(m::_HorizonRootModel,left,delta)
    target=_hr_basis(m,left).I0+delta
    # From a turning point sqrt(h) is analytic in phase, not a cancelled quadratic residual.
    w = left in m.turns ? abs(m.D * delta * _horizon_root_Sh(m.h2 * delta^2)) / 2 : nothing
    return _hr_basis_z(m,_hr_z_of_target(m,target);mino=target,sqrt_h=w)
end
# Mino time from the infinity endpoint (orbits that reach infinity, h0 >= 0)
function _hr_mino(m::_HorizonRootModel,r)
    m.h0>=0 || error("The horizon-root orbit does not reach infinity.")
    return -m.k0inf-_hr_basis(m,r).I0
end

function _radial_increment(model,left,right)
    o=zero(_float_type(left,right))
    left==right && return (I0=o,I1=o,I2=o,J1=o,J2=o)
    if right<left
        value=_radial_increment(model,right,left)
        return NamedTuple{keys(value)}(Tuple(-item for item in values(value)))
    end
    l=_hr_basis(model,left); r=_hr_basis(model,right)
    return (I0=r.I0-l.I0,I1=r.I1-l.I1,I2=r.I2-l.I2,
        J1=r.J1-l.J1,J2=r.J2-l.J2)
end

# Increment from `left` to the radius reached after the Mino-time step `delta`. Horizon-root
# models evaluate the endpoint in z = 1/(r - 1) directly (see _horizon_root_model); the others
# go through the radius.
function _coordinate_increment_from(model,a,energy,lz,q,left,delta,radius)
    l=_hr_basis(model,left); r=_hr_basis_from(model,left,delta)
    increment=(I0=r.I0-l.I0,I1=r.I1-l.I1,I2=r.I2-l.I2,J1=r.J1-l.J1,J2=r.J2-l.J2)
    return (mino=increment.I0,
        t=energy*increment.I2+2energy*increment.I1+
          3energy*increment.I0+4energy*increment.J1,
        phi=2a*energy*increment.J1,tau=increment.I2)
end

function _coordinate_increment(model,a,energy,lz,q,left,right)
    increment=_radial_increment(model,left,right)
    return (mino=increment.I0,
        t=energy*increment.I2+2energy*increment.I1+
          3energy*increment.I0+4energy*increment.J1,
        phi=2a*energy*increment.J1,tau=increment.I2)
end

_rstar(r)=r+2log((r-1)/2)-2/(r-1)
_phi_h(a,r)=-a/(r-1)

function _extremal_polar_engine(a,energy,lz,q,sector,phase;
        axis=nothing,polar_hemisphere::Symbol=:north)
    if axis!==nothing
        T=_float_type(a,energy,lz,q)
        o=zero(T)
        z0=axis===:north ? one(T) : axis===:south ? -one(T) :
            throw(ArgumentError("axis must be :north or :south"))
        formula=lambda -> (z=z0,uz=o,sin2=o,theta=acos(z0),phi=o,t=o,tau=lambda)
        return (formula=formula,sector=:axis_constant,phase=o,
            metadata=(sector=:axis_constant,phase=o,phase_convention=:not_applicable))
    end
    selected=sector!==nothing ? sector :
        _axis_crossing(a,energy,lz,q) ? :axis_crossing :
        iszero(q) ? :equatorial : q>0 ? :pendular : :vortical
    T=_float_type(a,energy,lz,q)
    solution=_polar_solution(a,energy,lz,q,selected,T(phase);
        hemisphere=polar_hemisphere)
    return (formula=solution.formula,sector=selected,phase=T(phase),metadata=solution.metadata,
        solution=solution)
end

function _check_domain(lambda,domain)
    low,high=domain.mino
    left_ok=domain.endpoint_closed[1] ? lambda>=low : lambda>low
    right_ok=domain.endpoint_closed[2] ? lambda<=high : lambda<high
    left_ok && right_ok || throw(DomainError(lambda,
        "Mino time lies outside $(domain.mino) with closure $(domain.endpoint_closed)."))
    return float(lambda)
end

function _component_kind(case_id)
    case_id===:A1 && return :periodic
    case_id in (:A2,:K1,:K3,:K6,:K9,:A_X2) && return :constant
    case_id===:K4 && return :homoclinic_outer
    case_id in (:B1,:B2,:B3,:B4,:B5,:B6,:B7,:B8,:B9) && return :horizon_to_turn
    case_id in (:K2,:K5,:K8,:K11) && return :horizon_to_repeated
    case_id in (:C1,:C2,:C3,:C4,:C5,:C6,:C7,:C8,:C9,:C10,:C11,:C12) &&
        return :direct_capture
    case_id in (:K7,:K10) && return :infinity_to_repeated
    case_id in (:D1,:D2) && return :scatter
    case_id in TRAPPED_CASE_IDS && return :trapped
    case_id in (:B_X1,:B_X2) && return :horizon_root_turn
    case_id in (:C_X2,:C_X4,
        :C_X1,:C_X3) && return :horizon_root_from_infinity
    case_id in (:D_X2,:D_X1) && return :horizon_root_scatter
    case_id===:A_X1 && return :horizon_root_island
    error("No exact-extremal motion kind is registered for $(case_id).")
end

function _strict_specs(a,energy,lz,q)
    kerr_polar_admissibility(a,energy,lz,q).admissible ||
        return (structure=kerr_geo_root_structure(a,energy,lz,q),specs=NamedTuple[],
            excluded=(:polar_motion_inadmissible,))
    if energy<0
        structure=kerr_geo_root_structure(a,energy,lz,q)
        ph=2energy-a*lz
        ph>0 || return (structure=structure,specs=NamedTuple[],
            excluded=(:trapped_horizon_root_or_past_oriented,))
        disposition=_trapped_disposition(energy,structure)
        disposition!==nothing ||
            return (structure=structure,specs=NamedTuple[],
                excluded=(:trapped_root_topology_not_class_n,))
        turn=first(structure.exterior).radius
        spec=(id=_trapped_case_id(disposition),broad=:trapped,
            formula=kerr_geo_case(_trapped_case_id(disposition)).FormulaFamily,
            disposition=disposition,lower=1.0,upper=turn,
            kind=:trapped,component=nothing)
        return (structure=structure,specs=[spec],excluded=())
    end
    classification=kerr_geo_classify(a,energy,lz,q)
    structure=classification.Status.root_structure
    specs=NamedTuple[]
    for component in classification.Components
        component.CaseId===nothing && continue
        id=component.CaseId
        lower=component.LowerEndpoint.Radius
        upper=component.UpperEndpoint.Radius
        kind=_component_kind(id)
        push!(specs,(id=id,broad=component.BroadClass,
            formula=component.FormulaFamily,lower=lower,upper=upper,
            kind=kind,component=component))
    end
    return (structure=structure,specs=specs,
        excluded=classification.ExcludedCaseIds,classification=classification)
end

# positive roots z = r − 1 of a z² + b z + c (the exterior roots at |a| = 1); `linear` for E = 1
function _positive_roots(a,b,c; linear)
    T=_float_type(a,b,c)
    exterior=_root_atol(T)+_root_rtol(T)
    if linear
        iszero(b) && return T[]
        root=-c/b
        return root>exterior ? [root] : T[]
    end
    disc=b^2-4a*c
    disc>=0 || return T[]
    roots=sort([(-b-sqrt(disc))/(2a),(-b+sqrt(disc))/(2a)])
    return [root for root in roots if root>exterior]
end

function _horizon_root_specs(a,energy,lz,q)
    structure=kerr_geo_root_structure(a,energy,lz,q)
    _horizon_root(a,energy,lz) || error("Horizon-root exact-extremal classification requires P_H=0.")
    energy>0 || return (structure=structure,specs=NamedTuple[],
        excluded=(:horizon_root_nonpositive_energy,))
    q>=0 || return (structure=structure,specs=NamedTuple[],
        excluded=(:horizon_root_negative_Q_polar_inadmissible,))
    x=energy^2; aa=_e2m1(energy); bb=4x-2; cc=3x-1-q
    tol=_horizon_root_coefficient_tol(_float_type(a,energy,lz,q))
    if abs(bb)<=tol && abs(cc)<=tol
        return (structure=structure,specs=NamedTuple[],excluded=(:horizon_root_quadruple_forbidden,))
    end
    regime=kerr_energy_regime(energy)
    roots=_positive_roots(aa,bb,cc;linear=regime===:parabolic)
    specs=NamedTuple[]
    if regime===:hyperbolic
        if cc>tol
            push!(specs,(id=:C_X2,broad=:capture,formula=:EXT_H2,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        elseif abs(cc)<=tol
            push!(specs,(id=:C_X4,broad=:capture,formula=:EXT_H3,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        else
            turn=only(roots)+1
            push!(specs,(id=:D_X2,broad=:scatter,formula=:EXT_H2_OUTER,
                lower=turn,upper=Inf,kind=:horizon_root_scatter,component=nothing))
        end
    elseif regime===:parabolic
        if cc>tol
            push!(specs,(id=:C_X1,broad=:capture,formula=:EXT_H2,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        elseif abs(cc)<=tol
            push!(specs,(id=:C_X3,broad=:capture,formula=:EXT_H3,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        else
            turn=only(roots)+1
            push!(specs,(id=:D_X1,broad=:scatter,
                formula=:EXT_H2_OUTER,lower=turn,upper=Inf,
                kind=:horizon_root_scatter,component=nothing))
        end
    else
        if cc>tol
            turn=only(roots)+1
            push!(specs,(id=:B_X1,broad=:plunge,formula=:EXT_H2,
                lower=1.0,upper=turn,kind=:horizon_root_turn,component=nothing))
        elseif abs(cc)<=tol && bb>tol
            turn=only(roots)+1
            push!(specs,(id=:B_X2,broad=:plunge,formula=:EXT_H3,
                lower=1.0,upper=turn,kind=:horizon_root_turn,component=nothing))
        elseif bb>tol && length(roots)==2
            qstable=x^2/(1-x)
            if abs(q-qstable)<=_tol(_float_type(q),1e-10)*max(1.0,abs(qstable))
                radius=1+(2x-1)/(1-x)
                push!(specs,(id=:A_X2,broad=:stable,
                    formula=:EXT_CONSTANT,lower=radius,upper=radius,
                    kind=:constant,component=nothing))
            elseif q<qstable
                push!(specs,(id=:A_X1,broad=:stable,
                    formula=:EXT_H2_ISLAND,lower=roots[1]+1,upper=roots[2]+1,
                    kind=:horizon_root_island,component=nothing))
            end
        end
    end
    isempty(specs) && return (structure=structure,specs=specs,
        excluded=(:horizon_root_no_future_exterior_component,))
    return (structure=structure,specs=specs,excluded=())
end

function _classification(a,energy,lz,q)
    (a==1 || a==-1) || throw(DomainError(a,
        "Exact-extremal dispatch requires a=+1 or a=-1 exactly."))
    ph=2energy-a*lz
    result=_horizon_root(a,energy,lz) ? _horizon_root_specs(a,energy,lz,q) :
        _strict_specs(a,energy,lz,q)
    metric_limit=a==1 ? :extremal_plus : :extremal_minus
    return merge(result,(metric_limit=metric_limit,P_H=ph,
        case_ids=Tuple(spec.id for spec in result.specs)))
end

function _default_reference(lower,upper)
    if isfinite(lower) && isfinite(upper)
        return lower+(upper-lower)/2
    elseif isfinite(lower)
        return lower+max(1.0,0.5abs(lower))
    end
    error("Cannot choose a reference radius without one finite endpoint.")
end

# Constant radius (circular or spherical orbit at a repeated root).
function _extremal_constant_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    radius=spec.lower
    metadata=_constant_radial_metadata(radius)
    radial_r=lambda -> radius
    radial_sign=lambda -> 0.0
    p=energy*(radius^2+1)-a*lz; y=radius-1
    tr=(radius^2+1)*p/y^2
    pr=a*p/y^2-a*energy
    radial_history=lambda -> (t=tr*lambda,phi=pr*lambda,tau=radius^2*lambda)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=metadata.domain, reference=metadata.reference)
end

# Periodic motion between two turning points; λ = 0 at the inner one.
function _extremal_periodic_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    lower=spec.lower; upper=spec.upper
    duration=_radial_increment(model,lower,upper).I0
    period=2duration
    folded = function (lambda)
        n=round(Int,lambda/period); rem=lambda-n*period
        outward=rem>=0; phase=abs(rem)
        radius=_hr_inverse_from(model,lower,phase)
        return (cycle=n,outward=outward,radius=radius,phase=phase)
    end
    radial_r=lambda -> folded(float(lambda)).radius
    radial_sign=lambda -> folded(float(lambda)).outward ? 1.0 : -1.0
    half=_coordinate_increment(model,a,energy,lz,q,lower,upper)
    radial_history = function (lambda)
        folded_value=folded(float(lambda))
        part=_coordinate_increment_from(model,a,energy,lz,q,lower,
            folded_value.phase,folded_value.radius)
        values=folded_value.outward ? part :
            (t=-part.t,phi=-part.phi,tau=-part.tau,mino=-part.mino)
        return (t=2folded_value.cycle*half.t+values.t,
            phi=2folded_value.cycle*half.phi+values.phi,
            tau=2folded_value.cycle*half.tau+values.tau)
    end
    domain=(mino=(-Inf,Inf),endpoint_closed=(false,false),
        endpoint_roles=(:infinite_past_worldline,:infinite_future_worldline),
        radial_period=period)
    reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
        t_phi_zero_lambda=0.0,t_phi_zero_radius=lower,tau_zero_event=:finite_turning_point,
        lambda_regular=nothing)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference)
end

# Motion symmetric about a turning point: scatter, homoclinic or horizon-root
# asymptotics on both sides.
function _extremal_two_sided_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    turn=spec.kind===:homoclinic_outer ? spec.upper :
        spec.kind===:horizon_root_turn ? spec.upper : spec.lower
    duration=kind in (:scatter,:horizon_root_scatter) ?
        _hr_mino(model,turn) : Inf
    is_scatter=kind in (:scatter,:horizon_root_scatter)
    radial_state = function (lambda)
        lam=float(lambda)
        abs(lam)<duration || throw(DomainError(lambda,
            "Mino time must lie inside the two-sided branch domain."))
        radius=_hr_inverse_from(model,turn,(is_scatter ? 1 : -1)*abs(lam))
        orientation=is_scatter ?
            (lam<0 ? -1.0 : lam>0 ? 1.0 : 0.0) :
            (lam<0 ? 1.0 : lam>0 ? -1.0 : 0.0)
        return (radius=radius,sign=orientation)
    end
    radial_r=lambda -> radial_state(lambda).radius
    radial_sign=lambda -> radial_state(lambda).sign
    radial_history = function (lambda)
        radial_value=radial_state(lambda)
        forward=_coordinate_increment_from(model,a,energy,lz,q,turn,
            (is_scatter ? 1 : -1)*abs(float(lambda)),radial_value.radius)
        increment=(t=-forward.t,phi=-forward.phi,tau=-forward.tau)
        orientation=is_scatter ?
            (lambda<0 ? 1.0 : lambda>0 ? -1.0 : 0.0) :
            (lambda<0 ? -1.0 : lambda>0 ? 1.0 : 0.0)
        return (t=orientation*increment.t,
            phi=orientation*increment.phi,
            tau=orientation*increment.tau)
    end
    domain=(mino=(-duration,duration),endpoint_closed=(false,false),
        endpoint_roles=kind in (:scatter,:horizon_root_scatter) ?
            (:past_infinity,:future_infinity) :
            kind===:horizon_root_turn ?
            (:past_horizon_root_asymptote,
             :future_horizon_root_asymptote) :
            (:past_repeated_root_asymptote,:future_repeated_root_asymptote))
    reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
        t_phi_zero_lambda=0.0,t_phi_zero_radius=turn,tau_zero_event=:finite_turning_point,
        lambda_regular=nothing)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference)
end

# From infinity to a repeated root (or a horizon root); λ = 0 at a reference radius.
function _extremal_infall_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    lower=spec.lower
    ref=reference_radius===nothing ? _default_reference(lower,spec.upper) :
        float(reference_radius)
    ref>lower || error("Reference radius must lie above the lower endpoint.")
    infinity_delta=_hr_mino(model,ref)
    lambda_min=-infinity_delta
    radial_state = function (lambda)
        lam=float(lambda); lam>lambda_min || throw(DomainError(lambda,
            "Mino time must exceed the past-infinity endpoint."))
        radius=_hr_inverse_from(model,ref,-lam)
        return (radius=radius,sign=-1.0)
    end
    radial_r=lambda -> radial_state(lambda).radius
    radial_sign=lambda -> -1.0
    radial_history = function (lambda)
        radius_value=radial_r(lambda)
        increment=_coordinate_increment_from(model,a,energy,lz,q,ref,-float(lambda),
            radius_value)
        return (t=-increment.t,phi=-increment.phi,tau=-increment.tau)
    end
    domain=(mino=(lambda_min,Inf),endpoint_closed=(false,false),
        endpoint_roles=(:past_infinity,
            kind===:horizon_root_from_infinity ? :future_horizon_root_asymptote :
            :future_repeated_root_asymptote))
    reference=(lambda0_event=:reference_radius,t_phi_zero_event=:reference_radius,
        t_phi_zero_lambda=0.0,t_phi_zero_radius=ref,tau_zero_event=:reference_radius,
        lambda_regular=nothing)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference)
end

# Leading Laurent/Puiseux behaviour for repeated (double or triple) horizon roots.
function _extremal_endpoint_metadata(kind,endpoint,model,a,energy)
    return if model isa _HorizonRootModel &&
            kind in (:horizon_root_turn,:horizon_root_from_infinity)
        c=model.h2; b=model.h1
        if abs(c)>_horizon_root_coefficient_tol(_float_type(c))
            (kind=:exact_extremal_horizon_root_double,
             horizon_multiplicity=2,series_type=:laurent,
             finite_regular_endpoint=false,
             v_kernel_leading=(exponent=-2,
                coefficient=4energy/sqrt(c)-2),
             psi_kernel_leading=(exponent=-2,
                coefficient=a*(2energy/sqrt(c)-1)),
             matching_rule=:horizon_root_asymptotic_no_finite_endpoint)
        else
            (kind=:exact_extremal_horizon_root_triple,
             horizon_multiplicity=3,series_type=:positive_y_puiseux,
             finite_regular_endpoint=false,
             v_kernel_leading=(exponent=-2.5,
                coefficient=4energy/sqrt(b)),
             psi_kernel_leading=(exponent=-2.5,
                coefficient=2a*energy/sqrt(b)),
             v_primitive_leading=(exponent=-1.5,
                coefficient=-8energy/(3sqrt(b))),
             psi_primitive_leading=(exponent=-1.5,
                coefficient=-4a*energy/(3sqrt(b))),
             matching_rule=:horizon_root_asymptotic_no_finite_endpoint)
        end
    else
        (kind=:no_horizon_endpoint,series_order=0,
         matching_rule=:not_applicable)
    end

end

function _make_trajectory(a,energy,lz,q,spec,structure;
        polar_sector=nothing,polar_phase=0.0,polar_hemisphere::Symbol=:north,
        axis=nothing,reference_radius=nothing)
    kind=spec.kind
    if kerr_geo_tier(spec.id)!==:extremal && kind!==:scatter
        return _extremal_engine_member(a,energy,lz,q,spec,structure;
            polar_sector=polar_sector,polar_phase=polar_phase,
            polar_hemisphere=polar_hemisphere,axis=axis,
            reference_radius=reference_radius)
    end
    # D1 and D2 never reach the horizon: their coordinates come from the radial engine, as at
    # |a| < 1 (the engine keeps t, φ, τ at rounding up to E → 1⁺ and far from the hole)
    if kind===:scatter
        axis===nothing || error("Axis trajectories use the axis-infall constructors; there is no scatter member on the axis.")
        polar=_extremal_polar_engine(a,energy,lz,q,polar_sector,polar_phase;
            polar_hemisphere=polar_hemisphere)
        return _scatter_member(a,energy,lz,q,spec.id,polar.solution,
            Tuple(item.radius for item in structure.real_roots);component=spec.component,
            structure=structure,
            formula_family=spec.formula,tier=:extremal,
            status=(metric_limit=a==1 ? :extremal_plus : :extremal_minus,
                metric_dispatch=a==1 ? :exact_positive : :exact_negative,motion_kind=kind,
                endpoint=_extremal_endpoint_metadata(kind,nothing,nothing,a,energy),
                polar_sector=polar.sector),
            chart=(rstar=_rstar,phi_h=r -> _phi_h(a,r)))
    end
    model=kind===:constant ? nothing : _horizon_root_model(energy,q;
        turns=Tuple(x for x in (spec.lower,spec.upper) if isfinite(x) && x!=1))
    pol=_extremal_polar_engine(a,energy,lz,q,polar_sector,polar_phase;
        axis=axis,polar_hemisphere=polar_hemisphere)
    polar=hasproperty(pol,:solution) ? pol.solution : pol
    builder=kind===:constant ? _extremal_constant_branch :
        kind===:horizon_root_island ? _extremal_periodic_branch :
        kind in (:horizon_root_scatter,:horizon_root_turn) ? _extremal_two_sided_branch :
        _extremal_infall_branch
    branch=builder(kind,spec,model,a,energy,lz,q,pol,reference_radius)
    check=lambda->_check_domain(lambda,branch.domain)
    radial=branch.radial_history
    p=_polar_primitive(polar)
    coords=(t=lambda->radial(lambda).t+p(lambda)[1],
        phi=lambda->radial(lambda).phi+p(lambda)[2],
        tau=lambda->radial(lambda).tau+p(lambda)[3])
    radial_r=branch.radial_r
    track=(r=lambda->radial_r(check(lambda)),check=check,check_bl=check,
        coords=coords,sign_r=branch.radial_sign,domain=branch.domain,
        reference=merge(branch.reference,(polar_phase=pol.metadata.phase,
            polar_phase_convention=pol.metadata.phase_convention)),trajectory=(;))
    return _engine_member(spec.broad,spec.id,a,energy,lz,q,polar,track,spec.component,structure,
        nothing,:extremal,
        (radial=structure.real_roots,raw=structure.raw_roots,lower=spec.lower,upper=spec.upper),
        (metric_limit=a==1 ? :extremal_plus : :extremal_minus,
            metric_dispatch=a==1 ? :exact_positive : :exact_negative,
            formula_family=spec.formula,motion_kind=kind,
            endpoint=_extremal_endpoint_metadata(kind,nothing,model,a,energy),
            disposition_id=get(spec,:disposition,nothing),polar_sector=pol.sector),
        (rstar=_rstar,phi_h=r->_phi_h(a,r)))
end

"""
    kerr_geo_extremal_family(a, E, Lz, Q; polar_sector=nothing, polar_phase=0.0,
                             polar_hemisphere=:north, axis=nothing, reference_radius=nothing)

Classify `(E, Lz, Q)` at exact `a = +1` or `a = -1` and construct every admitted member.
With `P_H = 2E − aLz = 0` the members are the horizon-root cases A-X1, A-X2, B-X1, B-X2,
C-X1…C-X4, D-X1, D-X2; otherwise they keep their primary case IDs. `polar_sector` selects the
polar sector, `polar_phase` is the polar phase at `λ = 0`, and `polar_hemisphere`
selects `:north` or `:south` for motion confined to one hemisphere.
`axis = :north` or `:south` puts
the motion on the spin axis (`Lz = 0`, `Q = a²(1 − E²)`), and `reference_radius` is the
radius at `λ = 0` for members that come in from infinity to a repeated or horizon root
(K7, K10, C-X1…C-X4).
"""
function kerr_geo_extremal_family(a::Real,energy::Real,lz::Real,q::Real;
        polar_sector=nothing,polar_phase=0.0,polar_hemisphere::Symbol=:north,
        axis=nothing,reference_radius=nothing)
    (a==1 || a==-1) || throw(DomainError(a,
        "Exact-extremal family dispatch requires exact a=+1 or a=-1."))
    if axis!==nothing
        _on_axis(a,energy,lz,q) || error("Axis initial data requires Lz = 0 and Q = a^2(1 - E^2).")
    end
    classification=_classification(a,float(energy),float(lz),float(q))
    members=Tuple(_make_trajectory(float(a),float(energy),float(lz),float(q),
        spec,classification.structure;polar_sector=polar_sector,
        polar_phase=polar_phase,polar_hemisphere=polar_hemisphere,
        axis=axis,reference_radius=reference_radius)
        for spec in classification.specs)
    metric_limit=a==1 ? :extremal_plus : :extremal_minus
    metric_dispatch=a==1 ? :exact_positive : :exact_negative
    status=(supported=!isempty(members),
        metric_dispatch=metric_dispatch,
        case_ids=classification.case_ids,excluded=classification.excluded)
    return KerrGeoExtremalFamily(metric_limit,
        (a=float(a),E=float(energy),Lz=float(lz),Q=float(q)),
        classification,members,status)
end

"""
    kerr_geo_extremal(a, E, Lz, Q; case_id=nothing, kwargs...)

Return the member `case_id` of `kerr_geo_extremal_family(a, E, Lz, Q; kwargs...)`. Without
`case_id` the constants must admit exactly one member.
"""
function kerr_geo_extremal(a::Real,energy::Real,lz::Real,q::Real;
        case_id=nothing,kwargs...)
    family=kerr_geo_extremal_family(a,energy,lz,q;kwargs...)
    if case_id===nothing
        length(family.Members)==1 || error(
            "Exact-extremal constants admit $(family.Classification.case_ids); provide case_id or use kerr_geo_extremal_family.")
        return only(family.Members)
    end
    members=filter(member -> member.CaseId===case_id,family.Members)
    length(members)==1 || error("Case $(case_id) is not an admitted exact-extremal member.")
    return only(members)
end
