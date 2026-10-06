# The radius and its horizon distance share one analytic radial phase. Keeping the
# distance separately is necessary when an exterior turning root rounds to r+.
struct _HorizonRelativeRadius{T,J}
    kind::Symbol
    parameters::NTuple{4,T}
    landen::J
    omega::T
    horizon_phase::T
    horizon_time::T
    horizon::Tuple{T,T}
    separation::T
    momentum::T
    turn_gap::T
    amplitude_denominator::T
    shifted::NTuple{5,Tuple{T,T}}
end

function _horizon_shifted_polynomial(a, E, L, Q)
    aa,ee,ll,qq = _wide.(float.((a,E,L,Q)))
    h,d = _wide_rplus(a)
    p = _wide_sub(_wide_mul(_wide(2.0),_wide_mul(ee,h)),_wide_mul(aa,ll))
    u = _wide_sub(ll,_wide_mul(aa,ee))
    hh = _wide_mul(h,h); e2 = _wide_mul(ee,ee); dd = _wide(d)
    k = _wide_add(_wide_add(hh,_wide_mul(u,u)),qq)
    times(n,x) = _wide_mul(_wide(float(n)),x)
    co = (_wide_mul(p,p),
        _wide_sub(times(4,_wide_mul(_wide_mul(ee,h),p)),_wide_mul(dd,k)),
        _wide_sub(_wide_sub(_wide_add(times(4,_wide_mul(e2,hh)),
            times(2,_wide_mul(ee,p))),k),times(2,_wide_mul(dd,h))),
        _wide_sub(_wide_sub(times(4,_wide_mul(e2,h)),times(2,h)),dd),
        _wide_sub(e2,_wide(1.0)))
    return h,d,p,co
end

function _horizon_turn_gap(co, horizon, turn)
    difference = _wide_sub(_wide(turn),horizon)
    delta = max(difference[1]+difference[2],zero(turn))
    derivative = _wide_derivative_coefficients(co,1)
    for _ in 1:64
        next = delta-_wide_evalpoly(delta,co)/_wide_evalpoly(delta,derivative)
        next == delta && break
        delta = next
    end
    delta > 0 || error("The simple radial turning point is not outside the horizon.")
    return delta
end

function _horizon_relative_model(a,E,L,Q,model)
    h,d,p,co = _horizon_shifted_polynomial(a,E,L,Q)
    delta = _horizon_turn_gap(co,h,model.turn)
    x = model.roots; kind = model.kind
    if kind in (:elliptic_double_d_below_simple,:elliptic_double_simple_below_d,
            :elliptic_triple_below_horizon)
        return _interior_repeated_horizon_model(E,model,h,d,p,co,delta)
    end
    distance(root)=begin
        value=_wide_sub(h,_wide(root))
        value[1]+value[2]
    end
    if kind === :four_simple_inner
        x1,x2,x3,x4 = x.x1,x.x2,x.x3,x.x4
        n = (x2-x1)/(x3-x1); n1 = (x3-x2)/(x3-x1)
        m = (x4-x3)*(x2-x1)/((x4-x2)*(x3-x1))
        m1 = (x3-x2)*(x4-x1)/((x4-x2)*(x3-x1))
        omega = sqrt(-_e2m1(E)*(x4-x2)*(x3-x1))/2
        denominator_h=n1*distance(x1)
        angle = atan(sqrt(delta),sqrt(denominator_h))
        parameters = (x3-x2,n,1.0,-n)
    elseif kind === :b2_outer_double_inner_simple
        span = x.x2-x.x1; total = x.x3-x.x1
        n1 = (x.x3-x.x2)/total
        omega = sqrt(-_e2m1(E)*total*(x.x3-x.x2))/2
        denominator_h=n1*distance(x.x1)
        angle = atan(sqrt(delta),sqrt(denominator_h))
        m = zero(span); m1 = one(span)
        parameters = (span,n1,one(span),n1-1)
    elseif kind === :four_simple_single_exterior
        x1,x2,x3,x4 = x.x1,x.x2,x.x3,x.x4
        h1 = (x4-x3)/(x3-x1)
        weight = (x4-x3)*(x4-x1)/(x3-x1)
        m = (x4-x3)*(x2-x1)/((x4-x2)*(x3-x1))
        m1 = (x3-x2)*(x4-x1)/((x4-x2)*(x3-x1))
        omega = sqrt(-_e2m1(E)*(x4-x2)*(x3-x1))/2
        denominator_h=distance(x3)*(x4-x1)/(x3-x1)
        angle = atan(sqrt(delta),sqrt(denominator_h))
        parameters = (weight,1.0,1.0,h1)
    elseif kind === :two_real_complex_pair
        aa,bb = x.A,x.B; span = x.x2-x.x1
        d1 = x.x1-x.rho; d2 = x.x2-x.rho
        wm1 = d1>0 ? x.eta^2/(bb+d1) : bb-d1
        wp1 = d1<0 ? x.eta^2/(bb-d1) : bb+d1
        wm2 = d2>0 ? x.eta^2/(aa+d2) : aa-d2
        wp2 = d2<0 ? x.eta^2/(aa-d2) : aa+d2
        m = (wm1-wm2)*(wp2-wp1)/(4aa*bb)
        m1 = (wm2+wp1)*(aa+bb+span)/(4aa*bb)
        omega = sqrt(-_e2m1(E)*aa*bb)
        denominator_h=distance(x.x1)
        angle = 2atan(sqrt(bb)*sqrt(delta),sqrt(aa)*sqrt(denominator_h))
        parameters = (aa,bb,span,0.0)
    elseif kind === :b5_parabolic_three_simple
        span = x.x2-x.x1; total = x.x3-x.x1
        m = span/total; m1 = (x.x3-x.x2)/total
        omega = sqrt(total/2)
        denominator_h=m1*distance(x.x1)
        angle = atan(sqrt(delta),sqrt(denominator_h))
        parameters = (span,m1,0.0,0.0)
    elseif kind === :b6_hyperbolic_four_simple
        x1,x2,x3,x4 = x.x1,x.x2,x.x3,x.x4
        n = (x3-x2)/(x4-x2); n1 = (x4-x3)/(x4-x2)
        m = (x3-x2)*(x4-x1)/((x4-x2)*(x3-x1))
        m1 = (x2-x1)*(x4-x3)/((x4-x2)*(x3-x1))
        omega = sqrt(_e2m1(E)*(x4-x2)*(x3-x1))/2
        denominator_h=n1*distance(x2)
        angle = atan(sqrt(delta),sqrt(denominator_h))
        parameters = (x3-x2,n1,1.0,-n)
    else
        error("No simple-root horizon chart for $(kind).")
    end
    landen = _landen(m,m1)
    phase = _ellip_f(angle,m1)
    T = _float_type(a,E,L,Q)
    return _HorizonRelativeRadius{T,typeof(landen)}(kind,parameters,landen,omega,phase,
        phase/omega,h,d,p[1]+p[2],delta,denominator_h,co)
end

function _horizon_relative_state(r::_HorizonRelativeRadius{T},lambda) where {T}
    lambda == 0 && return (gap=r.turn_gap,velocity=zero(T),chart=r)
    du = r.omega*(lambda-r.horizon_time)
    mid = _ellipj_reduced(r.horizon_phase+du/2,r.landen)
    step = _ellipj_reduced(du/2,r.landen)
    ju = _ellipj_reduced(r.omega*lambda,r.landen)
    jh = _ellipj_reduced(r.horizon_phase,r.landen)
    den = step[2]^2+mid[3]^2*step[1]^2
    difference = 2step[1]*mid[2]*mid[3]/den
    if r.kind === :two_real_complex_pair
        aa,bb,span = r.parameters
        D(j) = bb*(j[2]>=0 ? 1+j[2] : j[1]^2/(1-j[2]))+
            aa*(j[2]<=0 ? 1-j[2] : j[1]^2/(1+j[2]))
        dc = -2mid[1]*step[1]*mid[3]*step[3]/den
        gap = 2aa*bb*span*dc/(D(ju)*D(jh))
        velocity = -2aa*bb*span*r.omega*ju[1]*ju[3]/D(ju)^2
    elseif r.kind === :b5_parabolic_three_simple
        span,m1 = r.parameters
        dsd = 2step[1]*mid[2]*step[3]/(den*ju[3]*jh[3])
        gap = -span*m1*dsd*(ju[1]/ju[3]+jh[1]/jh[3])
        velocity = -2r.omega*span*m1*ju[1]*ju[2]/ju[3]^3
    else
        weight,n,base,h1 = r.parameters
        Du = base*ju[2]^2+(base+h1)*ju[1]^2
        Dh = base*jh[2]^2+(base+h1)*jh[1]^2
        gap = -weight*n*difference*(ju[1]+jh[1])/(Du*Dh)
        velocity = -2r.omega*weight*n*ju[1]*ju[2]*ju[3]/Du^2
    end
    return (gap=gap,velocity=velocity,chart=r)
end

(r::_HorizonRelativeRadius)(lambda) =
    r.horizon[1]+(r.horizon[2]+_horizon_relative_state(r,lambda).gap)
_radial_state(r::_HorizonRelativeRadius,lambda) = _horizon_relative_state(r,lambda)

function _horizon_lambda_of_gap(r::_HorizonRelativeRadius{T},gap) where {T}
    gap==r.turn_gap && return zero(T)
    iszero(gap) && return r.horizon_time
    displacement=r.turn_gap-gap
    if r.kind===:two_real_complex_pair
        aa,bb,span=r.parameters
        angle=2atan(sqrt(bb*displacement),sqrt(aa*(r.amplitude_denominator+gap)))
    elseif r.kind===:b5_parabolic_three_simple
        span,m1=r.parameters
        angle=atan(sqrt(displacement),sqrt(r.amplitude_denominator+m1*gap))
    else
        weight,n,base,h1=r.parameters
        angle=atan(sqrt(base*displacement),sqrt(r.amplitude_denominator+(base+h1)*gap))
    end
    return _ellip_f(angle,r.landen.m1)/r.omega
end

struct _InteriorRepeatedHorizonRadius{T}
    kind::Symbol
    span::T
    ratio::T
    omega::T
    horizon_phase::T
    tangent_h::T
    horizon_time::T
    simple_gap::T
    repeated_gap::T
    root_separation::T
    horizon::Tuple{T,T}
    separation::T
    momentum::T
    turn_gap::T
    shifted::NTuple{5,Tuple{T,T}}
end

function _interior_repeated_horizon_model(E,model,h,d,p,co,delta)
    roots=model.roots
    repeated=roots.repeated
    triple=model.kind===:elliptic_triple_below_horizon
    simple=triple ? repeated : roots.simple
    span=roots.outer-simple
    distance=_wide_sub(h,_wide(simple))
    simple_gap=distance[1]+distance[2]
    repeated_distance=_wide_sub(h,_wide(repeated))
    repeated_gap=repeated_distance[1]+repeated_distance[2]
    y_h=sqrt(delta/simple_gap)
    T=typeof(y_h)
    if triple
        ratio=one(T)
        omega=sqrt(-_e2m1(E))*span/2
        phase=y_h
        tangent=y_h
        root_separation=zero(T)
    else
        A=roots.outer-repeated
        B=simple-repeated
        root_separation=abs(B)
        ratio=A/abs(B)
        omega=sqrt(-_e2m1(E)*A*abs(B))/2
        tangent=y_h/sqrt(ratio)
        phase=B>0 ? atan(tangent) :
            asinh(sqrt(-B)*sqrt(delta)/(sqrt(span)*sqrt(repeated_gap)))
        tangent=B>0 ? tan(phase) : tanh(phase)
    end
    return _InteriorRepeatedHorizonRadius{T}(model.kind,span,ratio,omega,phase,
        tangent,phase/omega,simple_gap,repeated_gap,root_separation,
        h,d,p[1]+p[2],delta,co)
end

function _radial_state(r::_InteriorRepeatedHorizonRadius{T},lambda) where {T}
    lambda==0 && return (gap=r.turn_gap,velocity=zero(T),chart=r)
    u=r.omega*lambda
    offset=r.omega*(r.horizon_time-lambda)
    if r.kind===:elliptic_triple_below_horizon
        tangent=u
        difference=offset
        derivative=one(T)
    elseif r.kind===:elliptic_double_d_below_simple
        tangent=tan(u)
        step=tan(offset)
        difference=(1+tangent^2)*step/(1-tangent*step)
        derivative=1+tangent^2
    else
        tangent=tanh(u)
        step=tanh(offset)
        difference=sech(u)^2*step/(1+tangent*step)
        derivative=sech(u)^2
    end
    denominator=1+r.ratio*tangent^2
    horizon_denominator=1+r.ratio*r.tangent_h^2
    gap=r.span*r.ratio*difference*(r.tangent_h+tangent)/
        (denominator*horizon_denominator)
    velocity=-2r.span*r.ratio*tangent*derivative*r.omega/denominator^2
    return (;gap,velocity,chart=r)
end

(r::_InteriorRepeatedHorizonRadius)(lambda)=
    r.horizon[1]+(r.horizon[2]+_radial_state(r,lambda).gap)

function _horizon_lambda_of_gap(r::_InteriorRepeatedHorizonRadius{T},gap) where {T}
    iszero(gap) && return r.horizon_time
    gap==r.turn_gap && return zero(T)
    displacement=r.turn_gap-gap
    tangent=sqrt(displacement/(r.simple_gap+gap)/r.ratio)
    phase=r.kind===:elliptic_triple_below_horizon ? tangent :
        r.kind===:elliptic_double_d_below_simple ? atan(tangent) :
        asinh(sqrt(r.root_separation)*sqrt(displacement)/
            (sqrt(r.span)*sqrt(r.repeated_gap+gap)))
    return phase/r.omega
end

function _relative_rates(c,state,kind,hs)
    delta = state.gap; h = state.chart
    radius = h.horizon[1]+(h.horizon[2]+delta)
    P = h.momentum+c.E*delta*(2(h.horizon[1]+h.horizon[2])+delta)
    if kind === :horizon
        K = radius^2+(c.L-c.a*c.E)^2+c.Q
        D = P+hs*abs(state.velocity)
        return ((radius^2+c.a^2)*K/D,c.a*K/D-c.a*c.E,radius^2)
    end
    delta_r = delta*(delta+h.separation)
    return ((radius^2+c.a^2)*P/delta_r,c.a*P/delta_r-c.a*c.E,radius^2)
end

function _relative_logs(state)
    gap = state.gap; d = state.chart.separation
    gap > 0 || return (-oftype(gap, Inf),-oftype(gap, Inf))
    ratio_log = gap<d ? log(gap)-log(d)-log1p(gap/d) : -log1p(d/gap)
    return log((gap+d)/2),ratio_log
end
function _relative_rstar(state)
    h = state.chart; gap = state.gap
    radius = h.horizon[1]+(h.horizon[2]+gap)
    if iszero(h.separation)
        return radius+2log(gap)-2/gap-2log(oftype(gap, 2))
    end
    log_inner,log_ratio = _relative_logs(state)
    return radius+2log_inner+2(h.horizon[1]+h.horizon[2])/h.separation*log_ratio
end
function _relative_azimuth(a,state)
    iszero(a) && return zero(state.gap)
    iszero(state.chart.separation) && return -a/state.gap
    return a/state.chart.separation*_relative_logs(state)[2]
end

struct _ParabolicCriticalHorizonRadius{T}
    horizon::Tuple{T,T}
    separation::T
    momentum::T
    span::T
    simple::T
    repeated::T
    frequency::T
    angle::T
    tangent_h::T
    sech_h_squared::T
    shifted::NTuple{5,Tuple{T,T}}
end

function _parabolic_critical_horizon_model(a,E,L,Q,model)
    h,d,p,co=_horizon_shifted_polynomial(a,E,L,Q)
    s,rc=model.roots.x1,model.roots.x2
    span=rc-s
    gap=_wide_sub(h,_wide(s)); distant=_wide_sub(_wide(rc),h)
    tangent=sqrt((gap[1]+gap[2])/span)
    return _ParabolicCriticalHorizonRadius{typeof(tangent)}(h,d,p[1]+p[2],span,s,rc,sqrt(span/2),
        atanh(tangent),tangent,(distant[1]+distant[2])/span,co)
end

function _radial_state(r::_ParabolicCriticalHorizonRadius,lambda)
    du=-r.frequency*lambda
    step=tanh(du)
    difference=step*r.sech_h_squared/(1+r.tangent_h*step)
    tangent=tanh(r.angle+du)
    gap=r.span*difference*(tangent+r.tangent_h)
    velocity=-2r.span*r.frequency*tangent*sech(r.angle+du)^2
    return (;gap,velocity,chart=r)
end

(r::_ParabolicCriticalHorizonRadius)(lambda)=
    r.horizon[1]+(r.horizon[2]+_radial_state(r,lambda).gap)

function _parabolic_critical_lambda_of_radius(r::_ParabolicCriticalHorizonRadius,radius)
    radius==r.repeated && return -oftype(r.span, Inf)
    tangent=sqrt((radius-r.simple)/r.span)
    gap=_wide_sub(_wide(float(radius)),r.horizon)
    difference=(gap[1]+gap[2])/(r.span*(tangent+r.tangent_h))
    sech_squared=(r.repeated-radius)/r.span
    denominator=(sech_squared+tangent^2*r.sech_h_squared)/(1+tangent*r.tangent_h)
    return -atanh(difference/denominator)/r.frequency
end
