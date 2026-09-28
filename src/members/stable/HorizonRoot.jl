# Horizon-root stable members (P(r+) = 0): A-H1 (stable island), A-H2 (stable spherical).

function _horizon_stable_input(a,energy,lz,q)
    metric=kerr_metric_limit(a)
    metric in (:subextremal,:near_extremal) || throw(DomainError(
        a,"A subextremal horizon-root stable orbit requires 0<|a|<1."))
    0<energy<1 || throw(DomainError(
        energy,"A horizon-root stable orbit requires 0<E<1."))
    polar=kerr_polar_admissibility(a,energy,lz,q)
    polar.admissible || throw(DomainError(
        (energy,lz,q),"The polar potential is inadmissible."))
    horizons=kerr_horizons(a)
    pplus=kerr_radial_momentum(a,energy,lz,horizons.rplus)
    abs(pplus)<=2e-10*max(1.0,abs(energy),abs(lz)) || throw(DomainError(
        pplus,"A horizon-root Stable member requires P(r₊) = 0."))
    structure=kerr_geo_root_structure(a,energy,lz,q)
    length(structure.horizon_coincident)==1 || error(
        "A horizon-root Stable member requires exactly one root of R on the outer horizon.")
    only(structure.horizon_coincident).multiplicity==1 || error(
        "The subextremal outer-horizon root must be simple.")
    return horizons,structure
end

# A-H2: on the exterior (stable) double root, like the Critical members on their root.
function _horizon_spherical_member(a,energy,lz,q,structure,polar)
    radius=only(structure.exterior).radius
    kerr_radial_derivatives(a,energy,lz,q,radius).R2<0 || error(
        "The horizon-root constant-radius member must be radially stable.")
    return _engine_member(:stable,:A_H2,a,energy,lz,q,polar,
        _on_root_track(a,energy,lz,polar,radius);
        roots=(radial=Tuple(root.radius for root in structure.real_roots),),
        status=(formula_family=:HC_CONSTANT,formula_kind=:constant_radius,
            stability=:stable))
end

# A-H1: libration between the two exterior roots, radial t, φ, τ from the closed-form moments
# (the outer-horizon pole has zero residue); λ = 0 at the inner turning point.
function _horizon_island_member(a,energy,lz,q,horizons,structure,polar)
    radii=Tuple(root.radius for root in structure.real_roots)
    length(radii)==4 || error("The horizon-root stable island requires four real roots.")
    lower=structure.exterior[1].radius
    upper=structure.exterior[2].radius
    model=_outer_four_real_model(energy,radii)
    base=model.basis(lower)
    half_mino=model.basis(upper).I0-base.I0
    rminus=horizons.rminus
    function radial_increment(radius)
        value=model.basis(radius)
        I0,I1,I2=value.I0-base.I0,value.I1-base.I1,value.I2-base.I2
        Jminus=model.pole(rminus,radius)-model.pole(rminus,lower)
        return (t=energy*I2+2energy*I1+energy*(a^2+4-2horizons.rplus)*I0+
                  4horizons.rminus*energy*Jminus,
            phi=2a*energy*Jminus,tau=I2)
    end
    half=radial_increment(upper)
    period=2half_mino
    function radial_state(lambda)
        lam=_finite_mino(lambda)
        cycle=floor(Int,lam/period)
        remainder=lam-cycle*period
        outward=remainder<=half_mino
        radius=model.inverse_from(lower,outward ? remainder : period-remainder)
        part=radial_increment(radius)
        folded(h,x)=2cycle*h+(outward ? x : 2h-x)
        return (radius=radius,sign=outward ? 1.0 : -1.0,t=folded(half.t,part.t),
            phi=folded(half.phi,part.phi),tau=folded(half.tau,part.tau))
    end
    p(lambda)=polar.formula(lambda)
    coords=(t=λ->radial_state(λ).t+p(λ).t,phi=λ->radial_state(λ).phi+p(λ).phi,
        tau=λ->radial_state(λ).tau+p(λ).tau)
    track=(r=λ->radial_state(λ).radius,check=_finite_mino,check_bl=_finite_mino,
        coords=coords,sign_r=λ->radial_state(λ).sign,
        domain=(mino=(-Inf,Inf),endpoint_closed=(false,false),
            endpoint_roles=(:infinite_past_worldline,:infinite_future_worldline),
            radial_period=period),
        reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0,t_phi_zero_radius=lower,tau_zero_event=:finite_turning_point,
            lambda_regular=nothing),trajectory=(;))
    return _engine_member(:stable,:A_H1,a,energy,lz,q,polar,track;
        roots=(radial=radii,),
        status=(formula_family=:HC_FF01,formula_kind=:outer_four_real_libration,
            radial_period=period))
end

"""
    kerr_geo_horizon_stable(a, E, Lz, Q; polar_sector=nothing, polar_phase=0.0)

The Stable member of constants whose radial potential has a simple root on the outer horizon
(P(r₊) = 0; 0 < |a| < 1, 0 < E < 1): A-H1, the libration between the two exterior roots
(λ = 0 at the inner turning point), or A-H2, the stable spherical orbit on the exterior
double root. t, φ and τ vanish at λ = 0, where `polar_phase` is the polar phase.
"""
function kerr_geo_horizon_stable(a::Real,energy::Real,lz::Real,q::Real;
        polar_sector=nothing,polar_phase::Real=0.0)
    horizons,structure=_horizon_stable_input(a,energy,lz,q)
    sector=polar_sector===nothing ? _constants_polar_sector(energy,lz,q) : polar_sector
    polar=_polar_solution(a,energy,lz,q,sector,float(polar_phase))
    if length(structure.exterior)==2 &&
            all(root->root.multiplicity==1,structure.exterior)
        return _horizon_island_member(a,energy,lz,q,horizons,structure,polar)
    elseif length(structure.exterior)==1 &&
            only(structure.exterior).multiplicity==2
        return _horizon_spherical_member(a,energy,lz,q,structure,polar)
    end
    error("These constants have no horizon-root Stable member (A-H1 or A-H2).")
end
