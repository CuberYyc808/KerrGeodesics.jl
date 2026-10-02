# Primary cases at exact |a| = 1 use the shared radial models and coordinate
# engine. Only their reference events and chart conventions differ.
function _extremal_radial_model(a,E,L,Q,spec,structure;axis=nothing)
    id=spec.id
    if id===:A1
        r4,r3,r2,r1=Tuple(root.radius for root in structure.real_roots)
        return _libration_model(E,r1,r2,r3,r4)
    elseif id in (:B1,:B2,:B3,:B4,:B5,:B6)
        return _plunge_radial_model(id,E,structure)
    elseif id in TRAPPED_CASE_IDS
        return _trapped_radial_model(spec.disposition,E,structure)
    elseif id in (:K2,:K4,:K5,:K7,:K8,:K10,:K11)
        radii=id in (:K4,:K5) ? Tuple(_double_root_factorization(a,E,L,Q,
            id===:K4 ? spec.lower : spec.upper)) : Tuple(r.radius for r in structure.real_roots)
        return _critical_radial_model(id,E,radii)
    elseif id===:C1
        return _c1_radial_model(a,L,Q;structure=structure)
    elseif id===:C3
        return _c3_radial_model(a,E,L,Q;structure=structure)
    elseif id===:C5
        return _four_complex_model(_four_complex_parameters(E,structure))
    elseif id in (:C2,:C4)
        return _capture_analytic_model(id,E,Tuple(r.radius for r in structure.real_roots))
    end
    return interior_repeated_radial_model(id,a,E,L,Q,structure)
end

function _extremal_engine_member(a,E,L,Q,spec,structure;
        polar_sector=nothing,polar_phase=0.0,polar_hemisphere::Symbol=:north,
        axis=nothing,reference_radius=nothing)
    pol=_extremal_polar_engine(a,E,L,Q,polar_sector,polar_phase;
        axis=axis,polar_hemisphere=polar_hemisphere)
    polar=hasproperty(pol,:solution) ? pol.solution : pol
    potential=_radial_potential_from_roots(a,E,L,Q,structure)
    kind=spec.kind
    model=kind===:constant ? nothing :
        _extremal_radial_model(a,E,L,Q,spec,structure;axis=axis)
    track=if kind===:constant
        _on_root_track(a,E,L,polar,spec.lower)
    elseif kind===:periodic
        velocity=model.velocity
        coords=_engine_coordinates(a,E,L,Q,model.radius,_polar_primitive(polar);
            potential=potential,domain=(0.0,model.period),period=model.period,
            ends=(:turning,:turning),σ=1.0)
        (r=model.radius,check=_finite_mino,check_bl=_finite_mino,coords=coords,
            sign_r=lambda->sign(velocity(lambda)),
            domain=(mino=(-Inf,Inf),endpoint_closed=(false,false),
                endpoint_roles=(:infinite_past_worldline,:infinite_future_worldline),
                radial_period=model.period),
            reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
                t_phi_zero_lambda=0.0,t_phi_zero_radius=spec.lower,
                tau_zero_event=:finite_turning_point,lambda_regular=nothing),trajectory=(;))
    elseif kind===:homoclinic_outer
        _homoclinic_track(a,E,L,Q,model,polar,potential)
    elseif kind===:infinity_to_repeated
        ref=reference_radius===nothing ? _default_reference(spec.lower,spec.upper) : reference_radius
        tr=_infinity_track(a,E,L,Q,model,polar,1.0,ref,potential)
        merge(tr,(reference=merge(tr.reference,(lambda_regular=nothing,)),))
    elseif kind===:trapped
        _extremal_trapped_track(a,E,L,Q,spec,model,polar,potential)
    else
        _extremal_horizon_track(a,E,L,Q,spec,model,polar,potential)
    end
    track=merge(track,(reference=merge(track.reference,(polar_phase=pol.metadata.phase,
        polar_phase_convention=pol.metadata.phase_convention)),))
    endpoint=if kind in (:horizon_to_turn,:horizon_to_repeated,:direct_capture,:trapped)
        base=(kind=:exact_extremal_strict,coordinate_method=:radial_engine,
            finite_regular_endpoint=true)
        kind===:trapped ? merge(base,(
            horizon_chart=(past=:retarded_u_chi,future=:advanced_v_psi),
            endpoint_zero=(past=(u=0.0,chi=0.0),future=(v=0.0,psi=0.0)),
            bl_endpoint_behavior=(:t_diverges_at_both_horizons,:phi_diverges_at_both_horizons))) : base
    else
        _extremal_endpoint_metadata(kind,nothing,nothing,a,E)
    end
    return _engine_member(spec.broad,spec.id,a,E,L,Q,polar,track,spec.component,structure,
        potential,:extremal,
        (radial=structure.real_roots,raw=structure.raw_roots,lower=spec.lower,upper=spec.upper),
        (metric_limit=a==1 ? :extremal_plus : :extremal_minus,
            metric_dispatch=a==1 ? :exact_positive : :exact_negative,
            formula_family=spec.formula,motion_kind=kind,endpoint=endpoint,
            disposition_id=get(spec,:disposition,nothing),polar_sector=pol.sector),
        (rstar=_rstar,phi_h=r->_phi_h(a,r)))
end

function _extremal_horizon_track(a,E,L,Q,spec,model,polar,potential)
    relative=model.kind===:k8_parabolic_outer_double ?
        _parabolic_critical_horizon_model(a,E,L,Q,model) : nothing
    c=relative===nothing ? model.mino(1.0) : relative.angle/relative.frequency
    orientation=model.inward ? 1.0 : -1.0
    kind=spec.kind
    lo=kind===:horizon_to_repeated ? -Inf :
        kind===:direct_capture ? -c : -abs(model.mino(spec.upper)-c)
    domain=(mino=(lo,0.0),endpoint_closed=(kind===:horizon_to_turn,true),
        endpoint_roles=(kind===:horizon_to_repeated ? :past_repeated_root_asymptote :
            kind===:direct_capture ? :past_infinity : :finite_turning_point,:future_horizon),
        bl_mino=(lo,0.0))
    check=lambda->_check_domain(lambda,domain)
    model_radius,model_mino,upper=model.radius,model.mino,spec.upper
    function radius(lambda)
        lambda=check(lambda)
        iszero(lambda) && return 1.0
        lambda==lo && isfinite(upper) && return upper
        return relative===nothing ? model_radius(c+orientation*lambda) : relative(lambda)
    end
    lambda_of(r)=r==1.0 ? 0.0 : relative===nothing ? (model_mino(r)-c)/orientation :
        _parabolic_critical_lambda_of_radius(relative,r)
    ref=lambda_of(1.0+min(1.0,(spec.upper-1.0)/2))
    ends=(kind===:horizon_to_repeated ? :asymptote :
        kind===:direct_capture ? :infinity : :turning,:horizon)
    multiplicity=spec.component===nothing ? 2 : spec.component.UpperEndpoint.Multiplicity
    radial_track=relative===nothing ? radius : relative
    coords=_engine_coordinates(a,E,L,Q,radial_track,_polar_primitive(polar);
        potential=potential,domain=domain.mino,ends=ends,σ=-1.0,
        rd=spec.upper,multiplicity=multiplicity,λ_bl=ref,λ_tau=0.0,
        λ_regular=0.0,σ_regular=-1.0)
    # The strict convention fixes t and phi through the horizon-anchored ingoing
    # chart, rather than setting them to zero at an exterior event.
    rstar(lambda)=relative===nothing ? _rstar(radius(lambda)) :
        _relative_rstar(_radial_state(relative,lambda))
    phi_h(lambda)=relative===nothing ? _phi_h(a,radius(lambda)) :
        _relative_azimuth(a,_radial_state(relative,lambda))
    t(lambda)=iszero(lambda) ? Inf : _coords_v(coords,lambda)-rstar(lambda)
    phi(lambda)=iszero(lambda) ? a*Inf : _coords_psi(coords,lambda)-phi_h(lambda)
    exactcoords=merge(_coordinate_functions(coords),(t=t,phi=phi))
    extras=merge((lambda_of_radius=lambda_of,
        u=lambda->t(check(lambda))-rstar(lambda),
        chi=lambda->phi(check(lambda))-phi_h(lambda)),
        _radius_increments(coords,lambda_of,-1.0))
    relative===nothing || (extras=merge(extras,(rstar=lambda->rstar(check(lambda)),)))
    return (r=radius,check=check,check_bl=check,coords=exactcoords,sign_r=lambda->-1.0,
        radial_track=radial_track,
        domain=domain,
        reference=(lambda0_event=:future_horizon,t_phi_zero_event=:regular_chart_at_future_horizon,
            t_phi_zero_lambda=NaN,t_phi_zero_radius=NaN,tau_zero_event=:future_horizon,
            lambda_regular=0.0),trajectory=extras)
end

function _extremal_trapped_track(a,E,L,Q,spec,model,polar,potential)
    duration=model.mino(1.0)
    domain=(mino=(-duration,duration),endpoint_closed=(true,true),
        endpoint_roles=(:past_horizon,:future_horizon),bl_mino=(-duration,duration))
    check=lambda->_check_domain(lambda,domain)
    model_radius=model.radius
    radius(lambda)=abs(check(lambda))==duration ? 1.0 : model_radius(abs(lambda))
    coords=_engine_coordinates(a,E,L,Q,radius,_polar_primitive(polar);
        potential=potential,domain=domain.mino,ends=(:horizon,:horizon),turn=0.0,
        σ=-1.0,λ_bl=0.0,λ_regular=duration,σ_regular=-1.0)
    retarded=_regular_chart(coords,1.0,-duration)
    t(lambda)=abs(lambda)==duration ? sign(lambda)*Inf : _coords_t(coords,lambda)
    phi(lambda)=abs(lambda)==duration ? sign(lambda)*a*Inf : _coords_phi(coords,lambda)
    return (r=radius,check=check,check_bl=check,coords=merge(_coordinate_functions(coords),(t=t,phi=phi)),
        sign_r=lambda->-sign(check(lambda)),domain=domain,
        reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
            t_phi_zero_lambda=0.0,t_phi_zero_radius=spec.upper,tau_zero_event=:finite_turning_point,
            lambda_regular=duration),
        trajectory=(u=lambda->retarded(check(lambda))[1],chi=lambda->retarded(check(lambda))[2]))
end
