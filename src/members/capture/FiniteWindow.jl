# Finite-window capture API `kerr_geo_capture` (E >= 1, infinity -> future horizon) and the
# closed-form r(λ), λ(r) of the E = 1 (C1) and E > 1 (C3) captures it uses.

"""
    KerrGeoCapture

E ≥ 1 orbit falling from infinity into the future horizon. `Formula` is
`:hyperbolic_capture` (E > 1) or `:parabolic_capture` (E = 1). Mino time λ is zero on
the horizon; the trajectory gives r, θ, the regular horizon coordinates v and ψ, and
finite-window increments of the Boyer–Lindquist t and φ (which diverge on the horizon).
"""
struct KerrGeoCapture
    Formula::Symbol
    EnergyRegime::Symbol
    Outcome::Symbol
    OrbitalParameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    Parametrization::String
    Roots::NamedTuple
    ReferenceZero::NamedTuple
    Trajectory::Any
    Velocity::Any
    Potentials::NamedTuple
    Residuals::NamedTuple
    Status::NamedTuple
end

function _capture_domain(lambda_infinity)
    return (mino=(lambda_infinity, 0.0), endpoint_roles=(:past_infinity, :future_horizon),
            endpoint_closed=(false, true))
end

function _c3_complex_parameters(a, energy, lz, q; atol=1e-10, structure=nothing)
    roots_all = structure===nothing ? radial_roots_for_constants(a, energy, lz, q) : structure.raw_roots
    real_roots = Float64[]
    complex_roots = ComplexF64[]
    if structure===nothing
        for root in roots_all
            if abs(imag(root)) <= atol
                push!(real_roots, real(root))
            else
                push!(complex_roots, root)
            end
        end
    else
        append!(real_roots,(root.radius for root in structure.real_roots))
        append!(complex_roots,_nonreal_roots(structure))
    end
    length(real_roots) == 2 || return nothing
    length(complex_roots) == 2 || return nothing
    # (polished: near E = 1 the companion roots lose digits to the far root r1 ≈ −2/(E² − 1))
    coefficients = kerr_radial_coefficients(a, energy, lz, q)
    real_roots = structure===nothing ? sort!([_polish_root(coefficients, x) for x in real_roots]) : sort!(real_roots)
    upper = complex_roots[argmax(imag.(complex_roots))]
    structure===nothing && (upper=_polish_root(coefficients,upper))
    eta = abs(imag(upper))
    eta > 0 || return nothing
    rplus = _rplus(a)
    real_roots[1] < real_roots[2] < rplus || return nothing
    lead = _e2m1(energy)
    lead > 0 || return nothing
    return (
        energy=float(energy),
        lead=lead,
        r1=real_roots[1],
        r2=real_roots[2],
        rho=real(upper),
        eta=eta,
        rplus=rplus,
        roots_all=Tuple(roots_all),
        real_roots=Tuple(real_roots),
    )
end

function _c3_wide_shape(c)
    square(x)=_wide_mul(x,x)
    eta2=square(_wide(c.eta))
    d2=_wide_sub(_wide(c.r2),_wide(c.rho))
    d1=_wide_sub(_wide(c.r1),_wide(c.rho))
    b=_wide_sqrt(_wide_add(square(d2),eta2))
    d=_wide_sqrt(_wide_add(square(d1),eta2))
    numerator=_wide_sub(_wide_mul(d2,_wide_neg(d1)),eta2)
    product=_wide_mul(b,d)
    n=_wide_div(numerator,product)
    # (BC)^2 - numerator^2 = eta^2 (r2-r1)^2 retains the small modulus complement.
    magnitude=numerator[1]<0 ? _wide_neg(numerator) : numerator
    small=_wide_div(_wide_mul(eta2,square(_wide_sub(_wide(c.r2),_wide(c.r1)))),
        _wide_mul(_wide(2.0),_wide_mul(product,_wide_add(product,magnitude))))
    m,m1=numerator[1]<0 ?
        (_wide_mul(_wide(.5),_wide_sub(_wide(1.0),n)),small) :
        (small,_wide_mul(_wide(.5),_wide_add(_wide(1.0),n)))
    return (B=b,C=d,n=n,m=m,m1=m1)
end
_c3_shape(c)=map(x->x[1]+x[2],_c3_wide_shape(c))
_c3_scale(c,shape)=begin
    E=_wide(c.energy)
    lead=_wide_mul(_wide_sub(E,_wide(1.0)),_wide_add(E,_wide(1.0)))
    scale=_wide_sqrt(_wide_mul(lead,_wide_mul(shape.B,shape.C)))
    scale
end

# λ(r) of C3 measured from infinity (λ(∞) = 0): minus the Mino time from r to infinity, in
# Carlson's form with the complex pair ρ ± iη (`_mino_to_infinity`). At the horizon,
# real Legendre amplitudes avoid cancellation inside the complex Carlson arguments.
function _c3_lambda_of_r(c, r)
    r == c.rplus || return -_mino_to_infinity(c.lead,
        (c.r1, c.r2, complex(c.rho, c.eta), complex(c.rho, -c.eta)), r)
    return _c3_horizon_lambda(c,_wide(r))
end

# Jacobi addition evaluates the phase difference from the horizon directly.
function _c3_horizon_lambda(c,horizon,radius=Inf)
    shape = _c3_wide_shape(c)
    B,C,m,m1 = shape.B,shape.C,shape.m,shape.m1
    square(x)=_wide_mul(x,x)
    times(n,x)=_wide_mul(_wide(float(n)),x)
    value(x)=x[1]+x[2]
    h1=_wide_sub(horizon,_wide(c.r1)); h2=_wide_sub(horizon,_wide(c.r2))
    ratio=_wide_div(h2,h1)
    br=_wide_mul(C,ratio); den=_wide_add(B,br)
    bc=_wide_mul(B,C)
    sh=_wide_div(times(2,_wide_sqrt(_wide_mul(bc,ratio))),den)
    ch=_wide_div(_wide_sub(B,br),den)
    source_ratio,ratio_difference=if isinf(radius)
        (_wide(1.0),_wide_div(_wide_sub(_wide(c.r2),_wide(c.r1)),h1))
    else
        x=_wide(float(radius)); x1=_wide_sub(x,_wide(c.r1))
        (_wide_div(_wide_sub(x,_wide(c.r2)),x1),
            _wide_div(_wide_mul(_wide_sub(_wide(c.r2),_wide(c.r1)),
                _wide_sub(x,horizon)),_wide_mul(x1,h1)))
    end
    source_br=_wide_mul(C,source_ratio); source_den=_wide_add(B,source_br)
    si=_wide_div(times(2,_wide_sqrt(_wide_mul(bc,source_ratio))),source_den)
    ci=_wide_div(_wide_sub(B,source_br),source_den)
    di=_wide_sqrt(_wide_add(m1,_wide_mul(m,square(ci))))
    dh=_wide_sqrt(_wide_add(m1,_wide_mul(m,square(ch))))
    ti=_wide_sqrt(_wide_div(C,B)); root_ratio=_wide_sqrt(ratio)
    source_root=_wide_sqrt(source_ratio)
    tangent=_wide_div(_wide_mul(ti,ratio_difference),
        _wide_mul(_wide_add(source_root,root_ratio),
            _wide_add(_wide(1.0),_wide_mul(square(ti),_wide_mul(source_root,root_ratio)))))
    sin_difference=_wide_div(times(2,tangent),_wide_add(_wide(1.0),square(tangent)))
    denominator=_wide_add(square(di),_wide_mul(m,_wide_mul(square(si),square(ch))))
    product=_wide_mul(m,_wide_mul(_wide_mul(si,sh),_wide_mul(ci,ch)))
    dd=_wide_mul(di,dh)
    dn_numerator=product[1]>=0 ? _wide_add(dd,product) :
        _wide_div(_wide_mul(denominator,
            _wide_add(m1,_wide_mul(m,_wide_mul(square(ci),square(ch))))),_wide_sub(dd,product))
    sn_difference=_wide_div(_wide_mul(sin_difference,_wide_add(denominator,dn_numerator)),
        _wide_mul(_wide_add(di,dh),denominator))
    cn_difference=_wide_div(_wide_add(_wide_mul(ci,ch),
        _wide_mul(_wide_mul(si,sh),dd)),denominator)
    finite=_ellip_f(value(sn_difference),abs(value(cn_difference)),value(m1))
    difference=cn_difference[1]>=0 ? finite : 2*_ellip_k(value(m1))-finite
    result=_wide_div(_wide(-difference),_c3_scale(c,shape))
    return value(result)
end

struct _C3HorizonRadius{P,J}
    params::P
    landen::J
    B::Float64; C::Float64; m::Float64; m1::Float64
    scale::Float64
    sh::Float64; ch::Float64; dh::Float64
    si::Float64; ci::Float64; di::Float64
    base::Float64
    lambda_infinity::Float64
    horizon::Tuple{Float64,Float64}
    separation::Float64
    momentum::Float64
    shifted::NTuple{5,Tuple{Float64,Float64}}
end

function _c3_horizon_track(a,E,L,Q,c)
    h,d,p,co=_horizon_shifted_polynomial(a,E,L,Q)
    shape=_c3_shape(c)
    B,C,m,m1=shape.B,shape.C,shape.m,shape.m1
    h1=_wide_sub(h,_wide(c.r1)); h2=_wide_sub(h,_wide(c.r2))
    ratio=(h2[1]+h2[2])/(h1[1]+h1[2]); den=B+C*ratio
    sh=2sqrt(B*C*ratio)/den; ch=(B-C*ratio)/den
    dh=sqrt(m1+m*ch^2)
    si=2sqrt(B*C)/(B+C); ci=(B-C)/(B+C); di=sqrt(m1+m*ci^2)
    base=2B*C*(c.r2-c.r1)/((h1[1]+h1[2])*den)
    scale=_c3_scale(c,_c3_wide_shape(c))
    return _C3HorizonRadius(c,_landen(m,m1),B,C,m,m1,scale[1]+scale[2],
        sh,ch,dh,si,ci,di,base,_c3_horizon_lambda(c,h),h,d,p[1]+p[2],co)
end

function _radial_state(r::_C3HorizonRadius,lambda)
    iszero(lambda) && return (gap=0.0,velocity=-r.momentum,chart=r)
    L=r.landen; m=r.m
    s,cn,dn=_ellipj_reduced(-lambda*r.scale/2,L)
    w=r.dh^2+m*r.sh^2*cn^2
    sp=(r.sh*cn*dn+s*r.ch*r.dh)/w
    dp=(r.dh*dn-m*r.sh*r.ch*s*cn)/w
    difference=-2sp*s*dp*dn/(dp^2+m*sp^2*cn^2)
    si,ci,di=_ellipj_reduced((lambda-r.lambda_infinity)*r.scale/2,L)
    wi=r.di^2+m*r.si^2*ci^2
    spi=(r.si*ci*di-si*r.ci*r.di)/wi
    dpi=(r.di*di+m*r.si*r.ci*si*ci)/wi
    infinity_difference=2spi*si*dpi*di/(dpi^2+m*spi^2*ci^2)
    den=(r.B+r.C)*infinity_difference
    weight=2r.B*r.C*(r.params.r1-r.params.r2)
    gap=weight*difference/(r.base*den)
    sj,cj,dj=_ellipj_reduced(-lambda*r.scale,L)
    wd=r.dh^2+m*r.sh^2*cj^2
    snu=(r.sh*cj*dj+sj*r.ch*r.dh)/wd
    dnu=(r.dh*dj-m*r.sh*r.ch*sj*cj)/wd
    velocity=weight*r.scale*snu*dnu/den^2
    return (gap=gap,velocity=velocity,chart=r)
end
(r::_C3HorizonRadius)(lambda)=
    r.horizon[1]+(r.horizon[2]+_radial_state(r,lambda).gap)

function _capture_positive_radial_interval(r_left, r_right)
    r_left == r_right && return (r_left, r_right, 1)
    r_left < r_right && return (r_left, r_right, 1)
    return (r_right, r_left, -1)
end



# r at Mino time δ after the infinity endpoint, as a closure over the constants. With
# u = F(φ|m) (φ = 0 at r2; u is √(lead B C) times the Mino time from r2), φ = am u and
# ratio = tan²(φ/2) B/C = (r − r2)/(r − r1),
# r = (r2 − ratio r1)/(1 − ratio). At infinity ratio = 1 (tan²(φ∞/2) = C/B), and with the
# offset η = u∞ − u = δ √(lead B C)
#     D = cn u − cn u∞ = 2 sn P sn(η/2) dn P dn(η/2)/(1 − m sn²P sn²(η/2)),  P = u∞ − η/2,
#     1 + cn u = 2B/(B + C) + D,  1 − ratio = (B + C)/C · D/(1 + cn u),
#     ratio = (B/C)(1 − cn u)/(1 + cn u),
# every quantity a sum of positive terms or a product, including E → 1⁺ where r1 → −∞ and
# u∞ → 2K. The Jacobi functions at u∞ are exact (sn φ∞ = 2√(BC)/(B + C), cn φ∞ =
# (B − C)/(B + C)); those of P follow from η/2 by the addition formulas.
function _c3_radius_from_infinity(c)
    shape = _c3_shape(c)
    m, B, C = shape.m, shape.B, shape.C
    L = _landen(m, shape.m1)
    k = sqrt(c.lead * B * C)
    s∞ = 2 * sqrt(B * C) / (B + C)
    c∞ = (B - C) / (B + C)
    d∞ = sqrt(shape.m1 + m * c∞^2)
    return function (δ)
        s, cq, dq = _ellipj_reduced(δ * k / 2, L)
        # 1 − m s∞² s² = d∞² + m s∞² cn²u and 1 − m sn²P s² = dn²P + m sn²P cn²u: sums of positive
        # terms (as differences they lose a digit next to the horizon, where m → 1)
        w = d∞^2 + m * s∞^2 * cq^2
        snP = (s∞ * cq * dq - s * c∞ * d∞) / w                  # P = u∞ − η/2
        dnP = (d∞ * dq + m * s∞ * c∞ * s * cq) / w
        D = 2 * snP * s * dnP * dq / (dnP^2 + m * snP^2 * cq^2)  # cn u − cn u∞
        onep = 2B / (B + C) + D                                 # 1 + cn u
        gap = (B + C) / C * D / onep                            # 1 − ratio
        ratio = B / C * (2 - onep) / onep
        return (c.r2 - ratio * c.r1) / gap
    end
end

# The numerator is a Jacobi difference from the horizon; the denominator is a
# difference from infinity. Both small endpoint differences are retained directly.
function _c3_radius_from_horizon(c)
    shape = _c3_shape(c)
    B, C, m, m1 = shape.B, shape.C, shape.m, shape.m1
    L = _landen(m, m1)
    scale = sqrt(c.lead * B * C)
    ratio = (c.rplus - c.r2) / (c.rplus - c.r1)
    denominator = B + C * ratio
    sh = 2 * sqrt(B * C * ratio) / denominator
    ch = (B - C * ratio) / denominator
    dh = sqrt(m1 + m * ch^2)
    gap = (c.r2 - c.r1) / (c.rplus - c.r1)
    base = 2 * B * C * gap / denominator
    si = 2sqrt(B*C)/(B+C)
    ci = (B-C)/(B+C)
    di = sqrt(m1+m*ci^2)
    lambda_infinity = _c3_lambda_of_r(c,c.rplus)
    return function (lambda)
        iszero(lambda) && return c.rplus
        s, cn, dn = _ellipj_reduced(-lambda * scale / 2, L)
        w = dh^2 + m * sh^2 * cn^2
        sp = (sh * cn * dn + s * ch * dh) / w
        dp = (dh * dn - m * sh * ch * s * cn) / w
        difference = -2 * sp * s * dp * dn / (dp^2 + m * sp^2 * cn^2)
        sinf,cinf,dinf = _ellipj_reduced((lambda-lambda_infinity)*scale/2,L)
        winf = di^2+m*si^2*cinf^2
        spinf = (si*cinf*dinf-sinf*ci*di)/winf
        dpinf = (di*dinf+m*si*ci*sinf*cinf)/winf
        infinity_difference = 2spinf*sinf*dpinf*dinf/(dpinf^2+m*spinf^2*cinf^2)
        return c.rplus + 2 * B * C * (c.r1 - c.r2) * difference /
            (base * (B + C) * infinity_difference)
    end
end

# The radial model of C3: `radius(δ)` at Mino time δ after the infinity endpoint and its
# inverse `mino(r)`, with the parameters and root data of the finite-window API
function _c3_radial_model(a, energy, lz, q; structure=nothing)
    c = _c3_complex_parameters(a, energy, lz, q;structure=structure)
    c === nothing && return nothing
    return (kind=:c3_hyperbolic_two_real_complex, params=c, rplus=c.rplus,
        radius=_c3_radius_from_infinity(c), horizon_radius=_c3_radius_from_horizon(c),
        mino=r -> -_c3_lambda_of_r(c, r), inward=true,
        roots=(real=c.real_roots, complex=(rho=c.rho, eta=c.eta)),
        rootdata=(radial=c.real_roots, complex=(rho=c.rho, eta=c.eta), lead=c.lead,
            shape=_c3_shape(c)))
end

function _c1_one_real_parameters(a, lz, q; atol=1e-10, structure=nothing)
    roots_all = structure===nothing ? radial_roots_for_constants(a, 1.0, lz, q) : structure.raw_roots
    real_roots = Float64[]
    complex_roots = ComplexF64[]
    if structure===nothing
        for root in roots_all
            if abs(imag(root)) <= atol
                push!(real_roots, real(root))
            else
                push!(complex_roots, root)
            end
        end
    else
        append!(real_roots,(root.radius for root in structure.real_roots))
        append!(complex_roots,_nonreal_roots(structure))
    end
    sort!(real_roots)
    length(real_roots) == 1 || return nothing
    length(complex_roots) == 2 || return nothing
    upper = complex_roots[argmax(imag.(complex_roots))]
    eta = abs(imag(upper))
    eta > 0 || return nothing
    rplus = _rplus(a)
    real_roots[1] < rplus || return nothing
    pplus = rplus^2 + a^2 - a * lz
    pplus > 0 || return nothing
    return (
        x0=real_roots[1],
        rho=real(upper),
        eta=eta,
        rplus=rplus,
        pplus=pplus,
        roots_all=Tuple(roots_all),
            real_roots=Tuple(real_roots),
    )
end

function _c1_shape(c)
    d = c.x0 - c.rho
    B = hypot(d, c.eta)
    # (B-d)(B+d)=eta^2 retains the smaller parameter next to a repeated root.
    m = d > 0 ? c.eta^2 / (2B * (B+d)) : (B-d) / (2B)
    m1 = d < 0 ? c.eta^2 / (2B * (B-d)) : (B+d) / (2B)
    return (d=d, B=B, m=m, m1=m1)
end

function _c1_lambda_of_r(c, r)
    shape = _c1_shape(c)
    r > c.x0 || error("C1 radius must lie above the real cubic root.")
    # Measure from infinity using the complementary amplitude, without 2K - F.
    phi = 2atan(sqrt(shape.B / (r - c.x0)))
    return -_ellip_f(phi, shape.m1) / sqrt(2 * shape.B)
end

# The radial model of C1: `radius(δ)` at Mino time δ after the infinity endpoint (am(2K − u) =
# π − am(u): the amplitude measured from infinity is small there) and its inverse `mino(r)`
function _c1_radial_model(a, lz, q; structure=nothing)
    c = _c1_one_real_parameters(a, lz, q;structure=structure)
    c === nothing && return nothing
    shape = _c1_shape(c)
    L = _landen(shape.m, shape.m1)
    # tan(am u / 2) = sn u / (1 + cn u)
    function radius(δ)
        sn, cn, _ = _ellipj_reduced(δ * sqrt(2 * shape.B), L)
        return c.x0 + shape.B * ((1 + cn) / sn)^2
    end
    return (kind=:c1_parabolic_one_real_complex, params=c, rplus=c.rplus,
        radius=radius, mino=r -> -_c1_lambda_of_r(c, r), inward=true,
        roots=(real=c.real_roots, complex=(rho=c.rho, eta=c.eta)),
        rootdata=(radial=c.real_roots, complex=(rho=c.rho, eta=c.eta), shape=shape))
end

function _unsupported_capture(parameters, constants, outcome, reason)
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    return KerrGeoCapture(outcome.formula, outcome.energy_regime, outcome.outcome, parameters,
        constants, "Mino", (radial=outcome.roots, root_class=outcome.root_class),
        (kind=:future_horizon_regular_v, lambda_horizon=0.0),
        nothing, nothing,
        (radial=r -> kerr_radial_potential(a, energy, lz, q, r),
         polar_z=z -> kerr_polar_z_potential(a, energy, lz, q, z)),
        (radial=nothing, polar_z=nothing),
        (supported=false, reason=reason))
end

"""
    kerr_geo_capture(a, p, e, x; input=:apex, kwargs...)
    kerr_geo_capture(a, E, Lz, Q; input=:constants, polar_phase=0.0)
    kerr_geo_capture(a, constants; kwargs...)

Finite-window capture orbit (E ≥ 1, from infinity into the future horizon). Inclined
orbits need constants input; `polar_phase` is the polar argument at the horizon.
Constants without a capture component return a record with `Status.supported == false`.
APEX input whose periapsis p/(1 + e) lies outside the horizon throws an `ArgumentError`: a
capture orbit has no periapsis there.
"""
function kerr_geo_capture(a::Real, b::Real, c::Real, d::Real; input::Symbol=:apex,
        constants=nothing, polar_phase=nothing)
    if input === :constants
        constants_tuple = (E=b, Lz=c, Q=d)
        parameters = (a=a, input=:constants)
    elseif input === :apex
        # p/(1+e) is a root of R; outside the horizon it is a turning point that bounds every
        # motion from infinity (and below the separatrix the APEX constants may not exist)
        b / (1 + c) > _rplus(a) && throw(ArgumentError(
            "kerr_geo_capture: the periapsis p/(1+e) = $(b / (1 + c)) lies outside the horizon " *
            "r₊ = $(_rplus(a)); a capture orbit has no periapsis there. " *
            "Pass the constants of motion with input=:constants."))
        if constants === nothing
            c_ = kerr_geo_constants_of_motion(a, b, c, d)
            constants_tuple = (E=c_["E"], Lz=c_["Lz"], Q=c_["Q"])
        else
            constants_tuple = constants
        end
        parameters = (a=a, p=b, e=c, x=d, input=:apex)
    else
        error("Unknown capture input type. Use :apex or :constants.")
    end
    outcome = _outcome_at_infinity(a, constants_tuple.E, constants_tuple.Lz, constants_tuple.Q)
    outcome.outcome === :capture || return _unsupported_capture(parameters, constants_tuple,
        outcome, "Outcome is $(outcome.outcome), not capture.")
    return _capture_finite_window(parameters, constants_tuple, outcome;
        polar_phase=polar_phase === nothing ? 0.0 : polar_phase)
end

function kerr_geo_capture(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_capture(a, constants.E, constants.Lz, constants.Q; input=:constants, kwargs...)
end

function kerr_geo_capture(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_capture(a, constants[1], constants[2], constants[3]; input=:constants, kwargs...)
end

function Base.show(io::IO, kg::KerrGeoCapture)
    print(io, "KerrGeoCapture(", kg.Formula, ", constants=")
    show(io, kg.ConstantsOfMotion)
    print(io, ", supported=", kg.Status.supported, ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoCapture)
    println(io, "KerrGeoCapture (", kg.EnergyRegime, ")")
    _show_summary_field(io, "Parameters", kg.OrbitalParameters)
    _show_summary_field(io, "Constants", kg.ConstantsOfMotion)
    _show_summary_field(io, "Formula", kg.Formula)
    println(io, "  Trajectory = (t(lambda), r(lambda), theta(lambda), phi(lambda))")
    _show_summary_status(io, kg.Status)
end

# Radial closed forms (r(λ), λ(r)) of the two finite-window capture formulas:
# C3 (E > 1, two real roots inside the horizon plus a complex pair) and
# C1 (E = 1, one real root plus a complex pair).
_capture_radial_model(formula, a, energy, lz, q) = formula === :hyperbolic_capture ?
    _c3_radial_model(a, energy, lz, q) : _c1_radial_model(a, lz, q)

"""
Finite-window capture orbit for E > 1 (two real roots inside the horizon plus a complex
pair) or E = 1 (one real root plus a complex pair): radial motion from infinity into the
future horizon, λ = 0 on the horizon where the regular v and ψ vanish. Equatorial orbits
and generic inclined ones with a fixed `polar_phase` are supported.
"""
function _capture_finite_window(parameters, constants, outcome; polar_phase=0.0)
    formula = outcome.formula
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    unsupported(msg) = _unsupported_capture(parameters, constants, outcome, msg)
    inclined = !iszero(q)
    inclined && parameters.input !== :constants && return unsupported(
        "Inclined capture needs constants input (the APEX polar-phase convention is not defined).")
    radial = _capture_radial_model(formula, a, energy, lz, q)
    radial === nothing && return unsupported(
        "The radial roots do not have the ordering required by the $formula formula.")
    polar = _window_polar_motion(a, energy, lz, q, polar_phase)

    # λ = 0 on the horizon; the infinity endpoint precedes it by the Mino time from infinity
    relative=formula===:hyperbolic_capture ? _c3_horizon_track(a,energy,lz,q,radial.params) : nothing
    lambda_infinity = relative===nothing ? -radial.mino(radial.rplus) : relative.lambda_infinity
    R(r) = kerr_radial_potential(a, energy, lz, q, r)
    Θ(z) = kerr_polar_z_potential(a, energy, lz, q, z)
    r_of_lambda=relative===nothing ? (λ->λ==0 ? radial.rplus : radial.radius(λ-lambda_infinity)) : relative
    rdot(λ)=relative===nothing ? -sqrt(max(R(r_of_lambda(λ)),0.0)) : _radial_state(relative,λ).velocity
    # Mino time of radius r (λ = 0 on the horizon) and the polar increments between radii.
    function λ_of(r)
        r >= radial.rplus || throw(DomainError(r, "Capture radius must lie outside the future horizon."))
        r==radial.rplus && return 0.0
        return relative===nothing ? lambda_infinity+radial.mino(r) :
            _c3_horizon_lambda(radial.params,relative.horizon,r)
    end
    polar_phi(r1, r2) = polar.phi(λ_of(r1)) - polar.phi(λ_of(r2))
    polar_t(r1, r2) = polar.t(λ_of(r1)) - polar.t(λ_of(r2))
    # τ, v, ψ and the radial increments: radial spectral engine (infinity → horizon) +
    # polar primitive; τ, v, ψ vanish on the horizon (λ = 0)
    coords = _engine_coordinates(a, energy, lz, q, r_of_lambda, polar.primitive;
        potential=_coefficient_potential(a, energy, lz, q),
        domain=(lambda_infinity, 0.0), ends=(:infinity, :horizon), σ=-1.0,
        λ_bl=0.5 * lambda_infinity, λ_tau=0.0, λ_regular=0.0, σ_regular=-1.0)
    inc = _radius_increments(coords, λ_of, -1.0)
    radial_time, radial_phi, proper, radial_v, radial_psi = inc.radial_time_increment,
        inc.radial_phi_increment, inc.radial_proper_increment, inc.radial_v_increment,
        inc.radial_psi_increment
    tau(λ) = _coords_tau(coords, λ)
    v(λ) = _coords_v(coords, λ)
    psi(λ) = _coords_psi(coords, λ)
    rstar(λ) = λ >= -MINO_ENDPOINT_TOL ? NaN : kerr_rstar(a, r_of_lambda(λ))
    utheta(λ) = inclined ? -polar.uz(λ) / sqrt(max(1 - polar.z(λ)^2, 0.0)) : 0.0
    phase = _polar_phase_metadata(polar, polar_phase, :future_horizon_regular_endpoint)

    return KerrGeoCapture(formula, outcome.energy_regime, :capture, parameters, constants, "Mino",
        merge(radial.rootdata, (root_class=outcome.root_class,
                                polar=inclined ? polar.rootdata : nothing)),
        (kind=:future_horizon_regular_v, lambda_horizon=0.0,
         lambda_infinity=lambda_infinity, v_horizon=0.0, psi_horizon=0.0,
         rplus=radial.rplus, polar_phase=inclined ? polar_phase : nothing),
        (t=nothing, r=r_of_lambda, theta=polar.theta, phi=nothing, q=polar.z, tau=tau,
         rstar=rstar, u=nothing, v=v, psi=psi,
         radial_phi_increment=radial_phi, radial_time_increment=radial_time,
         radial_proper_increment=proper,
         total_phi_increment=(r1, r2) -> radial_phi(r1, r2) + polar_phi(r1, r2),
         total_time_increment=(r1, r2) -> radial_time(r1, r2) + polar_t(r1, r2),
         radial_psi_increment=radial_psi, radial_v_increment=radial_v,
         total_psi_increment=(r1, r2) -> radial_psi(r1, r2) + polar_phi(r1, r2),
         total_v_increment=(r1, r2) -> radial_v(r1, r2) + polar_t(r1, r2)),
        (ut=nothing, ur=rdot, utheta=utheta, uphi=nothing),
        (radial=R, polar_z=Θ),
        (radial=λ -> rdot(λ)^2 - R(r_of_lambda(λ)),
         polar_z=λ -> polar.uz(λ)^2 - Θ(polar.z(λ))),
        (supported=true, reason="ok", domain=_capture_domain(lambda_infinity),
         polar_phase=phase))
end

# ---------------------------------------------------------------------------------------
# Radial integrals over a C1 / C3 radial window by adaptive Gauss-Kronrod in
# u = log(r − r₊) (which absorbs the 1/(r − r₊) horizon pole of the BL rates). Used only
# by the closed-form members built outside the spectral engine: the E = 1 axis infall and
# the exact-extremal (a = 1) capture bases.
# ---------------------------------------------------------------------------------------
function _capture_sqrt_radial(c, r)
    complex_part = (r - c.rho)^2 + c.eta^2
    polynomial = haskey(c, :lead) ?
        c.lead * (r - c.r1) * (r - c.r2) * complex_part :
        2 * (r - c.x0) * complex_part
    return sqrt(max(polynomial, 0.0))
end

function _capture_radial_quadrature(integrand, c, r_left, r_right)
    left, right, sign = _capture_positive_radial_interval(r_left, r_right)
    left == right && return 0.0
    left >= c.rplus || error("Capture radial increment requires a window outside the future horizon.")
    rp = c.rplus
    function g(u)
        y = exp(u)
        r = rp + y
        return integrand(r) * y / _capture_sqrt_radial(c, r)
    end
    # left == r+ is allowed for integrands regular at the horizon (u -> -Inf)
    lower = left == rp ? -Inf : log(left - rp)
    # maxevals caps pathological inputs; a well-posed window converges in a few hundred
    value, _ = quadgk(g, lower, log(right - rp); rtol=1.0e-13, atol=0.0, maxevals=20_000)
    return sign * value
end
