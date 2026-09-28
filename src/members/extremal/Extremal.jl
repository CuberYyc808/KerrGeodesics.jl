# Exact-extremal (|a| = 1) members: classification, radial models and assembly, for the X-tier
# members A-X1, A-X2, B-X1, B-X2, C-X1..C-X4, D-X1, D-X2 and the primary cases at |a| = 1.
# The radial parts of t, φ, τ come from each radial model's Mino-time bases I0, I1, I2 and
# horizon poles, matched to a regular series at a simple horizon; the polar parts come from
# the polar engine.
#
# Zero conventions (`ReferenceZero`) follow the other members except where the |a| = 1
# construction anchors differently: every member that ends on the future horizon (B, C, K2,
# K5, K8, K11 at |a| = 1) has λ = 0 and v = ψ = 0 on the future horizon, where t and φ
# diverge, so t and φ are fixed by that chart (`t_phi_zero_event =
# :regular_chart_at_future_horizon`) rather than at a finite reference radius; A1 at |a| = 1
# has λ = 0 and t = φ = 0 at the inner turning point instead of APEX initial phases; members
# that do not reach a horizon use v = t + r_*, ψ = φ + φ_H with no additive shift
# (`lambda_regular = nothing`).

const EXACT_ENDPOINT_RETAINED_ORDER = 14
const EXACT_ENDPOINT_GENERATED_ORDER = 16
const EXACT_ENDPOINT_LAST_TERM_TARGET = 1.0e-15
const EXACT_ENDPOINT_OMITTED_TARGET = 5.0e-17
const EXACT_ENDPOINT_OFFSETS =
    (0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001)

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

function Base.show(io::IO, ::MIME"text/plain", family::KerrGeoExtremalFamily)
    println(io, "KerrGeoExtremalFamily(")
    print(io, "    MetricLimit = "); show(io, family.MetricLimit); println(io, ",")
    print(io, "    ConstantsOfMotion = "); show(io, family.ConstantsOfMotion); println(io, ",")
    print(io, "    CaseIds = "); show(io, family.Classification.case_ids); println(io, ",")
    print(io, "    MemberCount = "); show(io, length(family.Members)); println(io, ",")
    print(io, "    Status = "); show(io, family.Status); println(io)
    print(io, ")")
end

function _basis_wrapper(name, basis, pole, inverse_i0; infinity_i0=nothing)
    inverse_from(left, delta) = inverse_i0(basis(left).I0 + delta)
    return (name=name, basis=basis, pole=pole, inverse_from=inverse_from,
        infinity_i0=infinity_i0, horizon_root=false)
end

function _outer_four_real_model(energy, ascending_roots)
    x1, x2, x3, x4 = ascending_roots
    r1, r2, r3, r4 = x4, x3, x2, x1
    kappa = -_e2m1(energy)
    m = (r1-r2)*(r3-r4)/((r1-r3)*(r2-r4))
    n = (r1-r2)/(r1-r3)
    d = r2-r3
    scale = 2/sqrt(kappa*(r1-r3)*(r2-r4))
    phi_of_r(r) = asin(sqrt(clamp(
        (r1-r3)*(r-r2)/((r1-r2)*(r-r3)), 0.0, 1.0)))
    function basis(r)
        phi=phi_of_r(r)
        f=Elliptic.F(phi,m)
        pin=_pi_real(n,phi,m)
        j2=_j2_legendre(n,m,phi)
        return (I0=scale*f,
            I1=scale*(r3*f+d*pin),
            I2=scale*(r3^2*f+2r3*d*pin+d^2*j2))
    end
    function pole(h,r)
        phi=phi_of_r(r)
        nh=n*(r3-h)/(r2-h)
        f=Elliptic.F(phi,m)
        return scale*(f/(r3-h)-
            d*_pi_real(nh,phi,m)/((r2-h)*(r3-h)))
    end
    function inverse_i0(target)
        sn=Elliptic.Jacobi.sn(clamp(target/scale,0.0,Elliptic.K(m)),m)
        s2=sn^2
        return (r2-r3*n*s2)/(1-n*s2)
    end
    return _basis_wrapper(:extremal_four_real_outer,basis,pole,inverse_i0)
end

function _inner_four_real_model(energy, ascending_roots)
    yroots=Tuple(-ascending_roots[index] for index in 4:-1:1)
    outer=_outer_four_real_model(energy,yroots)
    function basis(r)
        value=outer.basis(-r)
        return (I0=-value.I0,I1=value.I1,I2=-value.I2)
    end
    pole(h,r)=outer.pole(-h,-r)
    inverse_i0(target)=-outer.inverse_from(-ascending_roots[2],-target)
    return _basis_wrapper(:extremal_four_real_inner,basis,pole,inverse_i0)
end

function _finite_complex_model(energy, real_roots, complex_root)
    r2,r1=real_roots
    rho=real(complex_root); eta=abs(imag(complex_root))
    aa=hypot(r1-rho,eta); bb=hypot(r2-rho,eta)
    kappa=-_e2m1(energy)
    xi=sqrt(kappa*aa*bb)
    m=((r1-r2)^2-(aa-bb)^2)/(4aa*bb)
    fpar=4aa*bb/(aa-bb)^2
    function psi_of_r(r)
        yr=(bb*(r1-r)-aa*(r-r2))/(bb*(r1-r)+aa*(r-r2))
        return pi/2+asin(clamp(yr,-1.0,1.0))
    end
    function basis(r)
        psi=psi_of_r(r)
        f=Elliptic.F(psi,m); e=Elliptic.E(psi,m)
        pi_f=elliptic_pi(-1/fpar,psi,m)
        i0=f/xi
        i1=(aa*r2-bb*r1)*i0/(aa-bb)+
            (aa+bb)*(r1-r2)*pi_f/(2*(aa-bb)*xi)+
            atan((r1-r2)*sin(psi)/sqrt(4aa*bb*(1-m*sin(psi)^2)))/sqrt(kappa)
        combo=aa^2+2r2^2-bb^2-2r1^2
        boundary=sqrt(aa*bb/kappa)*((aa+bb)/(aa-bb)+cos(psi))*
            sin(psi)*sqrt(1-m*sin(psi)^2)/(fpar+sin(psi)^2)
        angle_y=2sin(psi)*sqrt(max(
            fpar*(1-m*sin(psi)^2)*(1+fpar*m),0.0))
        angle_x=fpar-(1+2fpar*m)*sin(psi)^2
        angle=atan(angle_y,angle_x)
        i2=(aa*r2^2-bb*r1^2)*i0/(aa-bb)+sqrt(aa*bb/kappa)*e-
            (aa+bb)*combo*pi_f/(4*(aa-bb)*xi)+boundary-
            combo*angle/(4*(r1-r2)*sqrt(kappa))
        return (I0=-i0,I1=-i1,I2=-i2)
    end
    function pole(h,r)
        psi=psi_of_r(r)
        i0=Elliptic.F(psi,m)/xi
        d=-sqrt(4aa*bb*(r1-h)*(h-r2))/
            (aa*(h-r2)+bb*(r1-h))
        pi_h=elliptic_pi(1/d^2,psi,m)
        first=(aa-bb)*i0/(aa*(r2-h)-bb*(r1-h))
        second=-(r1-r2)*(bb*(r1-h)-aa*(h-r2))*pi_h/
            (2xi*(r1-h)*(h-r2)*(bb*(r1-h)+aa*(h-r2)))
        root_term=sqrt((r1-r2)/(kappa*(r1-h)*(h-r2)*
            (aa^2*(h-r2)+bb^2*(r1-h)-(r1-r2)*(r1-h)*(h-r2))))
        w=sqrt(max(1-m*sin(psi)^2,0.0)); droot=sqrt(max(1-d^2*m,0.0))
        numerator=(d*droot+w*sin(psi))^2+m*(d^2-sin(psi)^2)^2
        denominator=(d*droot-w*sin(psi))^2+m*(d^2-sin(psi)^2)^2
        return -(first+second-root_term*log(numerator/denominator)/4)
    end
    function inverse_i0(target)
        u=-xi*target
        sn=Elliptic.Jacobi.sn(u,m); cn=Elliptic.Jacobi.cn(u,m)
        return (2aa*bb*(r1+r2)+(aa-bb)*(aa*r2-bb*r1)*sn^2+
            2aa*bb*(r1-r2)*cn)/(4aa*bb+(aa-bb)^2*sn^2)
    end
    return _basis_wrapper(:extremal_two_real_complex,basis,pole,inverse_i0)
end

function _c1_model(a,lz,q)
    c=_c1_one_real_parameters(a,lz,q)
    c===nothing && error("Exact-extremal C1 requires one real radial root inside the horizon, one complex root pair and P(r+) > 0.")
    shape=_c1_shape(c); reference=1.1
    lambda_reference=_c1_lambda_of_r(c,reference)
    # I1, I2 and the horizon pole by quadrature (see _capture_radial_quadrature): the
    # Legendre-form closed forms lose their branch for some windows.
    basis(r)=(I0=_c1_lambda_of_r(c,r)-lambda_reference,
        I1=_capture_radial_quadrature(x->x,c,reference,r),
        I2=_capture_radial_quadrature(x->x^2,c,reference,r))
    pole(h,r)=_capture_radial_quadrature(x->1/(x-h),c,reference,r)
    function inverse_i0(target)
        absolute=target+lambda_reference
        amplitude=Elliptic.Jacobi.am(absolute*sqrt(2shape.B),shape.m)
        return c.x0+shape.B*tan(amplitude/2)^2
    end
    infinity_i0=_c1_lambda_infinity(c)-lambda_reference
    return _basis_wrapper(:extremal_parabolic_one_real_complex,basis,pole,
        inverse_i0;infinity_i0=infinity_i0)
end

function _c3_model(a,energy,lz,q)
    c=_c3_complex_parameters(a,energy,lz,q)
    c===nothing && error("Exact-extremal C3 requires E > 1, two real radial roots inside the horizon and one complex root pair.")
    shape=_c3_shape(c); reference=1.1
    lambda_reference=_c3_lambda_of_r(c,reference)
    basis(r)=(I0=_c3_lambda_of_r(c,r)-lambda_reference,
        I1=_capture_radial_quadrature(x->x,c,reference,r),
        I2=_capture_radial_quadrature(x->x^2,c,reference,r))
    pole(h,r)=_capture_radial_quadrature(x->1/(x-h),c,reference,r)
    function inverse_i0(target)
        absolute=target+lambda_reference
        phi=Elliptic.Jacobi.am(absolute*sqrt(c.lead*shape.B*shape.C),shape.m)
        ratio=tan(phi/2)^2*shape.B/shape.C
        return (c.r2-ratio*c.r1)/(1-ratio)
    end
    infinity_i0=_c3_lambda_infinity(c)-lambda_reference
    return _basis_wrapper(:extremal_hyperbolic_two_real_complex,basis,pole,
        inverse_i0;infinity_i0=infinity_i0)
end

function _d1_model(radii)
    roots=(x1=radii[1],x2=radii[2],x3=radii[3])
    function basis(r)
        phi=_d1_amplitude(roots,r)
        return (I0=-_d1_i0_primitive(roots,phi),
            I1=-_d1_i1_primitive(roots,phi),
            I2=-_d1_i2_primitive(roots,phi))
    end
    pole(h,r)=-_d1_pole_primitive(
        roots,h,_d1_amplitude(roots,r))
    function inverse_i0(target)
        f=-target/_d1_scale(roots)
        sn=Elliptic.Jacobi.sn(f,_d1_modulus(roots))
        return roots.x1+_d1_A(roots)/sn^2
    end
    return _basis_wrapper(:extremal_parabolic_scatter,basis,pole,inverse_i0;
        infinity_i0=0.0)
end

function _d2_model(energy,radii)
    roots=(rA=radii[1],rB=radii[2],rC=radii[3],rD=radii[4])
    function basis(r)
        phi=_d2_amplitude(roots,r)
        return (I0=_d2_i0_primitive(energy,roots,phi),
            I1=_d2_i1_primitive(energy,roots,phi),
            I2=_d2_i2_primitive(energy,roots,phi))
    end
    pole(h,r)=_d2_pole_primitive(
        1.0,energy,roots,h,_d2_amplitude(roots,r))
    function inverse_i0(target)
        prefactor=_d2_prefactor(energy,roots)
        sn=Elliptic.Jacobi.sn(target/prefactor,_d2_modulus(roots))
        s2=sn^2; alpha=roots.rC-roots.rA; beta=roots.rD-roots.rA
        return (s2*beta*roots.rC-alpha*roots.rD)/(s2*beta-alpha)
    end
    infinity_i0=_d2_lambda_infinity(energy,roots)
    return _basis_wrapper(:extremal_hyperbolic_scatter,basis,pole,
        inverse_i0;infinity_i0=infinity_i0)
end

function _extremal_four_complex_model(a,energy,lz,q)
    parameters=_four_complex_parameters(energy,kerr_geo_root_structure(a,energy,lz,q))
    source=_four_complex_model(parameters)
    return _basis_wrapper(:extremal_four_complex,source.basis,source.pole,
        source.inverse_i0;infinity_i0=source.infinity_i0)
end

function _axis_extremal_model(energy)
    source=_axis_kerr_radial_model(1.0,energy)
    function horizon_pole(r)
        amplitude=source.chi(r)
        s=sin(amplitude)
        w=sqrt(max(1-source.m*s^2,0.0))
        logarithm=0.5log(abs((w-1)/(w+1)))
        f=_axis_legendre_f(amplitude,source.m)
        return (logarithm-f)/(2source.omega)
    end
    pole(h,r)=h==1.0 ? horizon_pole(r) : source.pole(h,r)
    return _basis_wrapper(:extremal_axis_legendre,source.basis,pole,
        source.radius_from_primitive;
        infinity_i0=source.lambda_primitive(Inf))
end

function _homoclinic_models(energy,radii)
    x1,rc,ra=radii
    source=_homoclinic_model(energy,x1,rc,ra)
    outer=_basis_wrapper(:extremal_homoclinic_outer,source.basis,source.pole,
        source.inverse_i0)
    inner=_basis_wrapper(:extremal_homoclinic_inner,source.inner_basis,
        source.inner_pole,source.inner_inverse_i0)
    return outer,inner
end

function _horizon_root_sqrt(value)
    value>=-1.0e-12 || throw(DomainError(value,
        "The horizon-root quadratic lies outside its allowed radial interval."))
    return sqrt(max(value,0.0))
end

# k0, k1, k2 are ∫ z^k dz/√h, h = h2 z² + h1 z + h0. At a turning point (root = true: z is a
# root of h) √h is exactly 0 and the inverse trigonometric function takes its branch value:
# from the rounded root they would carry an error ~ √eps (the square-root singularity).
_horizon_root_radicand(value, root) = root ? 0.0 : _horizon_root_sqrt(value)

function _horizon_root_k0(h2,h1,h0,z; root=false)
    d=h1^2-4h2*h0
    if h2>1e-14
        return root ? 0.0 : acosh(max((2h2*z+h1)/_horizon_root_sqrt(d),1.0))/sqrt(h2)
    elseif h2 < -1e-14
        arg=(2h2*z+h1)/_horizon_root_sqrt(d)
        return -asin(root ? sign(arg) : clamp(arg,-1.0,1.0))/sqrt(-h2)
    end
    return 2*_horizon_root_radicand(h1*z+h0,root)/h1
end

function _horizon_root_k1(h2,h1,h0,z; root=false)
    if abs(h2)>1e-14
        return _horizon_root_radicand(h2*z^2+h1*z+h0,root)/h2-
            h1*_horizon_root_k0(h2,h1,h0,z;root=root)/(2h2)
    end
    w=_horizon_root_radicand(h1*z+h0,root)
    return 2*(w^3/3-h0*w)/h1^2
end

function _horizon_root_k2(h2,h1,h0,z; root=false)
    if abs(h2)>1e-14
        h=h2*z^2+h1*z+h0
        k0=_horizon_root_k0(h2,h1,h0,z;root=root); k1=_horizon_root_k1(h2,h1,h0,z;root=root)
        s0=(2h2*z+h1)*_horizon_root_radicand(h,root)/(4h2)+
            (4h2*h0-h1^2)*k0/(8h2)
        return (s0-h1*k1-h0*k0)/h2
    end
    w=_horizon_root_radicand(h1*z+h0,root)
    return 2*(w^5/5-2h0*w^3/3+h0^2*w)/h1^3
end

function _horizon_root_model(energy,q; turns=())
    h0=_e2m1(energy); h1=4energy^2-2; h2=3energy^2-1-q
    # The natural variable is z = 1/(r - 1). Near the horizon root r = 1 + 1/z loses
    # all of z's digits (and saturates at r = 1), so the Mino-time path works in z:
    # basis_from(left, delta) evaluates the primitives at the z reached after delta.
    # (1/z is a root of the reversed quadratic exactly when z is a root of h)
    function basis_z(z; root=false)
        k0=_horizon_root_k0(h2,h1,h0,z;root=root)
        km1=-_horizon_root_k0(h0,h1,h2,inv(z);root=root)
        km2=-_horizon_root_k1(h0,h1,h2,inv(z);root=root)
        return (I0=-k0,I1=-(k0+km1),I2=-(k0+2km1+km2),
            J1=-_horizon_root_k1(h2,h1,h0,z;root=root),
            J2=-_horizon_root_k2(h2,h1,h0,z;root=root))
    end
    # the branch's turning points (the same radii its spec was built with) are roots of h
    basis(r)=basis_z(inv(r-1); root=r in turns)
    # z(k0) inverts k0(z); when the orbit reaches infinity (h0 >= 0) it is written as the
    # difference from z(k0∞) = 0, so z -> 0 (r -> infinity) keeps its digits
    d=h1^2-4h2*h0
    k0inf=h0>=0 ? _horizon_root_k0(h2,h1,h0,0.0) : NaN
    function z_of_target(target)
        k0=-target
        if h0>=0
            return if h2>1e-14
                s=sqrt(h2); sqrt(d)*sinh(s*(k0+k0inf)/2)*sinh(s*(k0-k0inf)/2)/h2
            elseif h2 < -1e-14
                s=sqrt(-h2); -sqrt(d)*cos(s*(k0+k0inf)/2)*sin(s*(k0-k0inf)/2)/h2
            else
                h1*(k0-k0inf)*(k0+k0inf)/4
            end
        end
        return if h2>1e-14
            (sqrt(d)*cosh(sqrt(h2)*k0)-h1)/(2h2)
        elseif h2 < -1e-14
            (-sqrt(d)*sin(sqrt(-h2)*k0)-h1)/(2h2)
        else
            ((h1*k0/2)^2-h0)/h1
        end
    end
    inverse_from(left,delta)=1+inv(z_of_target(basis(left).I0+delta))
    basis_from(left,delta)=basis_z(z_of_target(basis(left).I0+delta))
    infinity_i0=h0>=0 ? -_horizon_root_k0(h2,h1,h0,0.0) : nothing
    return (name=h2>1e-14 ? :extremal_horizon_root_cosh :
        (h2 < -1e-14 ? :extremal_horizon_root_cos : :extremal_horizon_root_linear),
        basis=basis,pole=nothing,inverse_from=inverse_from,basis_from=basis_from,
        infinity_i0=infinity_i0,horizon_root=true,h=(h2=h2,h1=h1,h0=h0))
end

function _strict_model(case_id,a,energy,lz,q,structure; disposition=nothing)
    radii=Tuple(root.radius for root in structure.real_roots)
    if case_id===:A1
        return _outer_four_real_model(energy,radii)
    elseif case_id===:B1
        return _inner_four_real_model(energy,radii)
    elseif case_id in (:B2,:K2,:K8,:B5,:K11,:B6)
        source=case_id in (:K2,:K8,:K11) ? _critical_radial_model(case_id,energy,radii) :
            _plunge_analytic_model(case_id,energy,radii)
        return _basis_wrapper(source.kind,source.basis,source.pole,
            source.inverse_i0;
            infinity_i0=haskey(source,:infinity_i0) ? source.infinity_i0 : nothing)
    elseif case_id===:B3
        return _outer_four_real_model(energy,radii)
    elseif case_id in TRAPPED_CASE_IDS
        disposition===:NFD01 && return _inner_four_real_model(energy,radii)
        disposition===:NFD02 && begin
            source=_plunge_analytic_model(:B2,energy,radii)
            return _basis_wrapper(source.kind,source.basis,source.pole,
                source.inverse_i0)
        end
        disposition===:NFD03 && return _outer_four_real_model(energy,radii)
        disposition===:NFD04 && begin
            upper=only(root for root in structure.raw_roots if imag(root)>0)
            return _finite_complex_model(energy,(radii[1],radii[2]),upper)
        end
        disposition===:NFD05 && begin
            source=_plunge_analytic_model(:B5,energy,radii)
            return _basis_wrapper(source.kind,source.basis,source.pole,
                source.inverse_i0)
        end
        disposition===:NFD06 && begin
            source=_plunge_analytic_model(:B6,energy,radii)
            return _basis_wrapper(source.kind,source.basis,source.pole,
                source.inverse_i0)
        end
        error("No exact-extremal Class N radial model is registered for disposition $(disposition).")
    elseif case_id===:B4
        upper=only(root for root in structure.raw_roots if imag(root)>0)
        return _finite_complex_model(energy,(radii[1],radii[2]),upper)
    elseif case_id in (:C2,:K7,:C4,:K10)
        source=case_id in (:K7,:K10) ? _critical_radial_model(case_id,energy,radii) :
            _capture_analytic_model(case_id,energy,radii)
        return _basis_wrapper(source.kind,source.basis,source.pole,
            source.inverse_i0;infinity_i0=source.infinity_i0)
    elseif case_id===:C1
        return _c1_model(a,lz,q)
    elseif case_id===:C3
        return _c3_model(a,energy,lz,q)
    elseif case_id===:C5
        is_axis=abs(lz)<=1.0e-12 &&
            abs(q-kerr_axis_carter_q(a,energy))<=
                1.0e-10*max(1.0,abs(q))
        return is_axis ? _axis_extremal_model(energy) :
            _extremal_four_complex_model(a,energy,lz,q)
    elseif case_id===:D1
        return _d1_model(radii)
    elseif case_id===:D2
        return _d2_model(energy,radii)
    elseif case_id in (:K4,:K5)
        outer,inner=_homoclinic_models(energy,radii)
        return case_id===:K4 ? outer : inner
    elseif case_id in (:B7,:B8,:B9,:C6,:C7,:C8,:C9,:C10,:C11,:C12)
        source=interior_repeated_radial_model(
            case_id,a,energy,lz,q,structure)
        return _basis_wrapper(source.kind,source.basis,source.pole,
            source.inverse_i0;
            infinity_i0=haskey(source,:infinity_i0) ? source.infinity_i0 : nothing)
    end
    error("No exact-extremal radial model is registered for $(case_id).")
end

function _radial_increment(model,left,right)
    left==right && return (I0=0.0,I1=0.0,I2=0.0,J1=0.0,J2=0.0)
    if right<left
        value=_radial_increment(model,right,left)
        return NamedTuple{keys(value)}(Tuple(-item for item in values(value)))
    end
    l=model.basis(left); r=model.basis(right)
    if model.horizon_root
        return (I0=r.I0-l.I0,I1=r.I1-l.I1,I2=r.I2-l.I2,
            J1=r.J1-l.J1,J2=r.J2-l.J2)
    end
    return (I0=r.I0-l.I0,I1=r.I1-l.I1,I2=r.I2-l.I2,
        J1=model.pole(1.0,right)-model.pole(1.0,left),J2=NaN)
end

function _strict_second_pole(a,energy,lz,q,left,right,increment)
    h=1.0; c3=2.0; c4=_e2m1(energy)
    polynomial=c4*increment.I2+c3*increment.I1/2-
        (c3*h/2+c4*h^2)*increment.I0
    rright=kerr_radial_potential(a,energy,lz,q,right)
    rleft=kerr_radial_potential(a,energy,lz,q,left)
    tolerance=1.0e-9*max(1.0,right^4)
    rright>=-tolerance || throw(DomainError(rright,
        "The right radial endpoint lies outside the allowed interval."))
    rleft>=-tolerance || throw(DomainError(rleft,
        "The left radial endpoint lies outside the allowed interval."))
    boundary=sqrt(max(rright,0.0))/(right-h)-
        sqrt(max(rleft,0.0))/(left-h)
    derivatives=kerr_radial_derivatives(a,energy,lz,q,h)
    return (polynomial-derivatives.R1*increment.J1/2-boundary)/derivatives.R
end

# Increment from `left` to the radius reached after the Mino-time step `delta`. Horizon-root
# models evaluate the endpoint in z = 1/(r - 1) directly (see _horizon_root_model); the others
# go through the radius.
function _coordinate_increment_from(model,a,energy,lz,q,left,delta,radius)
    model.horizon_root || return _coordinate_increment(model,a,energy,lz,q,left,radius)
    l=model.basis(left); r=model.basis_from(left,delta)
    increment=(I0=r.I0-l.I0,I1=r.I1-l.I1,I2=r.I2-l.I2,J1=r.J1-l.J1,J2=r.J2-l.J2)
    return (mino=increment.I0,
        t=energy*increment.I2+2energy*increment.I1+
          3energy*increment.I0+4energy*increment.J1,
        phi=2a*energy*increment.J1,tau=increment.I2)
end

function _coordinate_increment(model,a,energy,lz,q,left,right)
    increment=_radial_increment(model,left,right)
    j2=model.horizon_root ? increment.J2 :
        _strict_second_pole(a,energy,lz,q,left,right,increment)
    ph=2energy-a*lz
    if model.horizon_root
        return (mino=increment.I0,
            t=energy*increment.I2+2energy*increment.I1+
              3energy*increment.I0+4energy*increment.J1,
            phi=2a*energy*increment.J1,tau=increment.I2)
    end
    return (mino=increment.I0,
        t=energy*increment.I2+2energy*increment.I1+
          (5energy-a*lz)*increment.I0+(8energy-2a*lz)*increment.J1+2ph*j2,
        phi=2a*energy*increment.J1+a*ph*j2,tau=increment.I2)
end

_rstar(r)=r+2log((r-1)/2)-2/(r-1)
_phi_h(a,r)=-a/(r-1)

function _strict_endpoint_coefficients(a,energy,lz,q)
    order=EXACT_ENDPOINT_GENERATED_ORDER+2
    ph=2energy-a*lz
    ph>0 || error("Strict extremal endpoint coefficients require P_H>0.")
    p=zeros(Float64,order+1); p[1]=ph; p[2]=2energy; p[3]=energy
    sgeom=zeros(Float64,order+1)
    sgeom[1]=1+(energy-ph)^2+q; sgeom[2]=2; sgeom[3]=1
    rseries=_series_mul(p,p,order)
    for n in 2:order
        rseries[n+1]-=sgeom[n-1]
    end
    sqrt_r=zeros(Float64,order+1); sqrt_r[1]=ph
    for n in 1:order
        middle=sum(sqrt_r[k+1]*sqrt_r[n-k+1] for k in 1:(n-1);init=0.0)
        sqrt_r[n+1]=(rseries[n+1]-middle)/(2ph)
    end
    invsqrt=zeros(Float64,order+1); invsqrt[1]=inv(ph)
    for n in 1:order
        invsqrt[n+1]=-sum(sqrt_r[k+1]*invsqrt[n-k+1]
            for k in 1:n)/ph
    end
    w=_series_mul(p,invsqrt,order); w[1]=1.0; w[2]=0.0
    v=zeros(Float64,EXACT_ENDPOINT_GENERATED_ORDER+1)
    psi=similar(v)
    for n in 0:EXACT_ENDPOINT_GENERATED_ORDER
        wn=n==0 ? 0.0 : w[n+1]
        v[n+1]=2w[n+3]+2w[n+2]+wn
        psi[n+1]=a*(w[n+3]-energy*invsqrt[n+1])
    end
    return (v=v,psi=psi)
end

function _endpoint_match(coefficients)
    for y in EXACT_ENDPOINT_OFFSETS
        last=abs(coefficients[15]*y^15/15)
        omitted=abs(coefficients[16]*y^16/16)+
            abs(coefficients[17]*y^17/17)
        last<=EXACT_ENDPOINT_LAST_TERM_TARGET &&
            omitted<=EXACT_ENDPOINT_OMITTED_TARGET &&
            return (offset=y,last_retained=last,omitted_estimate=omitted)
    end
    error("Exact-extremal endpoint series did not meet its matching target.")
end

function _series_integral(coefficients,y)
    return sum(coefficients[n+1]*y^(n+1)/(n+1)
        for n in 0:EXACT_ENDPOINT_RETAINED_ORDER)
end

function _strict_endpoint_engine(model,a,energy,lz,q)
    coefficients=_strict_endpoint_coefficients(a,energy,lz,q)
    vm=_endpoint_match(coefficients.v); pm=_endpoint_match(coefficients.psi)
    function evaluate(component,r)
        r>=1 || throw(DomainError(r,"Exact-extremal endpoint radius must satisfy r>=1."))
        r==1 && return 0.0
        data=component===:v ? coefficients.v : coefficients.psi
        match=component===:v ? vm : pm
        y=r-1
        y<=match.offset && return _series_integral(data,y)
        rmatch=1+match.offset
        value=_series_integral(data,match.offset)
        increment=_coordinate_increment(model,a,energy,lz,q,rmatch,r)
        subtraction=component===:v ? _rstar(r)-_rstar(rmatch) :
            _phi_h(a,r)-_phi_h(a,rmatch)
        return value+(component===:v ? increment.t : increment.phi)-subtraction
    end
    return (v=r->evaluate(:v,r),psi=r->evaluate(:psi,r),
        order=EXACT_ENDPOINT_RETAINED_ORDER,generated_order=EXACT_ENDPOINT_GENERATED_ORDER,
        v_match=vm,psi_match=pm)
end

function _extremal_polar_engine(a,energy,lz,q,sector,phase;axis=nothing)
    if axis!==nothing
        z0=axis===:north ? 1.0 : axis===:south ? -1.0 :
            throw(ArgumentError("axis must be :north or :south"))
        formula=lambda -> (z=z0,uz=0.0,sin2=0.0,theta=acos(z0),phi=0.0,t=0.0,tau=lambda)
        return (formula=formula,sector=:axis_constant,phase=0.0,
            metadata=(sector=:axis_constant,phase=0.0,phase_convention=:not_applicable))
    end
    selected=sector!==nothing ? sector :
        abs(lz)<=1e-12 && q+a^2*_e2m1(energy)>1e-12 ? :axis_crossing :
        abs(q)<=1e-12 ? :equatorial : q>0 ? :pendular : :vortical
    solution=_polar_solution(a,energy,lz,q,selected,float(phase);hemisphere=:north)
    return (formula=solution.formula,sector=selected,phase=float(phase),metadata=solution.metadata)
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
    structure=kerr_geo_root_structure(a,energy,lz,q)
    kerr_polar_admissibility(a,energy,lz,q).admissible ||
        return (structure=structure,specs=NamedTuple[],
            excluded=(:polar_motion_inadmissible,))
    if energy<0
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

function _positive_roots(a,b,c)
    if abs(a)<=1e-14
        abs(b)>1e-14 || return Float64[]
        root=-c/b
        return root>0 ? [root] : Float64[]
    end
    disc=b^2-4a*c
    disc>=0 || return Float64[]
    roots=sort([(-b-sqrt(disc))/(2a),(-b+sqrt(disc))/(2a)])
    return [root for root in roots if root>1e-12]
end

function _horizon_root_specs(a,energy,lz,q)
    ph=2energy-a*lz
    abs(ph)<=1e-11 || error("Horizon-root exact-extremal classification requires P_H=0.")
    energy>0 || return (structure=nothing,specs=NamedTuple[],
        excluded=(:horizon_root_nonpositive_energy,))
    q>=0 || return (structure=nothing,specs=NamedTuple[],
        excluded=(:horizon_root_negative_Q_polar_inadmissible,))
    x=energy^2; aa=_e2m1(energy); bb=4x-2; cc=3x-1-q
    if abs(x-0.5)<=1e-11 && abs(q-0.5)<=1e-11
        return (structure=nothing,specs=NamedTuple[],excluded=(:horizon_root_quadruple_forbidden,))
    end
    roots=_positive_roots(aa,bb,cc)
    specs=NamedTuple[]
    if x>1+1e-12
        if cc>1e-11
            push!(specs,(id=:C_X2,broad=:capture,formula=:EXT_H2,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        elseif abs(cc)<=1e-11
            push!(specs,(id=:C_X4,broad=:capture,formula=:EXT_H3,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        else
            turn=only(roots)+1
            push!(specs,(id=:D_X2,broad=:scatter,formula=:EXT_H2_OUTER,
                lower=turn,upper=Inf,kind=:horizon_root_scatter,component=nothing))
        end
    elseif abs(x-1)<=1e-12
        if cc>1e-11
            push!(specs,(id=:C_X1,broad=:capture,formula=:EXT_H2,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        elseif abs(cc)<=1e-11
            push!(specs,(id=:C_X3,broad=:capture,formula=:EXT_H3,
                lower=1.0,upper=Inf,kind=:horizon_root_from_infinity,component=nothing))
        else
            turn=only(roots)+1
            push!(specs,(id=:D_X1,broad=:scatter,
                formula=:EXT_H2_OUTER,lower=turn,upper=Inf,
                kind=:horizon_root_scatter,component=nothing))
        end
    else
        if cc>1e-11
            turn=only(roots)+1
            push!(specs,(id=:B_X1,broad=:plunge,formula=:EXT_H2,
                lower=1.0,upper=turn,kind=:horizon_root_turn,component=nothing))
        elseif abs(cc)<=1e-11 && x>0.5+1e-12
            turn=only(roots)+1
            push!(specs,(id=:B_X2,broad=:plunge,formula=:EXT_H3,
                lower=1.0,upper=turn,kind=:horizon_root_turn,component=nothing))
        elseif x>0.5+1e-12 && length(roots)==2
            qstable=x^2/(1-x)
            if abs(q-qstable)<=1e-10*max(1.0,abs(qstable))
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
    isempty(specs) && return (structure=nothing,specs=specs,
        excluded=(:horizon_root_no_future_exterior_component,))
    return (structure=nothing,specs=specs,excluded=())
end

function _classification(a,energy,lz,q)
    (a==1 || a==-1) || throw(DomainError(a,
        "Exact-extremal dispatch requires a=+1 or a=-1 exactly."))
    ph=2energy-a*lz
    result=abs(ph)<=1e-11 ? _horizon_root_specs(a,energy,lz,q) :
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
    domain=(mino=(-Inf,Inf),endpoint_closed=(false,false),
        endpoint_roles=(:infinite_past_worldline,:infinite_future_worldline))
    radial_r=lambda -> radius
    radial_sign=lambda -> 0.0
    p=energy*(radius^2+1)-a*lz; y=radius-1
    tr=(radius^2+1)*p/y^2
    pr=a*p/y^2-a*energy
    radial_history=lambda -> (t=tr*lambda,phi=pr*lambda,tau=radius^2*lambda)
    reference=(lambda0_event=:polar_phase_reference,t_phi_zero_event=:polar_phase_reference,
        t_phi_zero_lambda=0.0,t_phi_zero_radius=radius,tau_zero_event=:polar_phase_reference,
        lambda_regular=nothing)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference)
end

# Periodic motion between two turning points; λ = 0 at the inner one.
function _extremal_periodic_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    lower=spec.lower; upper=spec.upper
    duration=_radial_increment(model,lower,upper).I0
    period=2duration
    folded = function (lambda)
        n=floor(Int,lambda/period); rem=lambda-n*period
        outward=rem<=duration; phase=outward ? rem : 2duration-rem
        radius=model.inverse_from(lower,phase)
        return (cycle=n,outward=outward,radius=radius,phase=phase)
    end
    radial_r=lambda -> folded(float(lambda)).radius
    radial_sign=lambda -> folded(float(lambda)).outward ? 1.0 : -1.0
    half=_coordinate_increment(model,a,energy,lz,q,lower,upper)
    radial_history = function (lambda)
        folded_value=folded(float(lambda))
        part=_coordinate_increment(model,a,energy,lz,q,lower,folded_value.radius)
        values=folded_value.outward ? part :
            (t=2half.t-part.t,phi=2half.phi-part.phi,
             tau=2half.tau-part.tau,mino=2half.mino-part.mino)
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
        model.infinity_i0-model.basis(turn).I0 : Inf
    is_scatter=kind in (:scatter,:horizon_root_scatter)
    radial_state = function (lambda)
        lam=float(lambda)
        abs(lam)<duration || throw(DomainError(lambda,
            "Mino time must lie inside the two-sided branch domain."))
        radius=model.inverse_from(turn,-abs(lam))
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
            -abs(float(lambda)),radial_value.radius)
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
    infinity_delta=model.infinity_i0-model.basis(ref).I0
    lambda_min=-infinity_delta
    radial_state = function (lambda)
        lam=float(lambda); lam>lambda_min || throw(DomainError(lambda,
            "Mino time must exceed the past-infinity endpoint."))
        radius=model.inverse_from(ref,-lam)
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

# Into the future horizon (from a turning point, a repeated root or infinity); λ = 0 on
# the horizon.
function _extremal_horizon_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    upper=spec.upper
    horizon_i0=model.basis(1.0).I0
    duration=kind===:horizon_to_repeated ? Inf :
        (isfinite(upper) ? model.basis(upper).I0-horizon_i0 :
         model.infinity_i0-horizon_i0)
    lambda_min=-duration
    radial_state = function (lambda)
        lam=float(lambda)
        lambda_min<lam<=0 || (lam==lambda_min && isfinite(upper)) ||
            throw(DomainError(lambda,"Mino time lies outside the member's domain, which ends on the future horizon at λ = 0."))
        radius=lam==0 ? 1.0 : model.inverse_from(1.0,-lam)
        return (radius=radius,sign=-1.0)
    end
    radial_r=lambda -> radial_state(lambda).radius
    radial_sign=lambda -> lambda==0 ? -1.0 : -1.0
    radial_history=lambda -> error("Horizon branches use regular endpoint coordinates.")
    past_role=kind===:horizon_to_repeated ? :past_repeated_root_asymptote :
        (isfinite(upper) ? :finite_turning_point : :past_infinity)
    domain=(mino=(lambda_min,0.0),
        endpoint_closed=(kind===:horizon_to_turn,true),
        endpoint_roles=(past_role,:future_horizon),
        bl_mino=(lambda_min,0.0))
    reference=(lambda0_event=:future_horizon,t_phi_zero_event=:regular_chart_at_future_horizon,
        t_phi_zero_lambda=NaN,t_phi_zero_radius=NaN,tau_zero_event=:future_horizon,
        lambda_regular=0.0)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference)
end

# Trapped between the past and future horizons, turning once in between.
function _extremal_trapped_branch(kind, spec, model, a, energy, lz, q, polar, reference_radius)
    turn=spec.upper
    duration=model.basis(turn).I0-model.basis(1.0).I0
    trapped_past_polar=polar.formula(-duration)
    trapped_future_polar=polar.formula(duration)
    trapped_radial_tau=model.basis(turn).I2-model.basis(1.0).I2
    radial_state = function (lambda)
        lam=float(lambda); -duration<=lam<=duration ||
            throw(DomainError(lambda,"Mino time lies outside the trapped domain between the past and future horizons."))
        delta=duration-abs(lam)
        radius=abs(lam)==duration ? 1.0 : model.inverse_from(1.0,delta)
        return (radius=radius,sign=lam<0 ? 1.0 : lam>0 ? -1.0 : 0.0)
    end
    radial_r=lambda -> radial_state(lambda).radius
    radial_sign=lambda -> radial_state(lambda).sign
    radial_history = function (lambda)
        lam=float(lambda); radial_value=radial_state(lam)
        orientation=Base.sign(lam)
        abs(lam)==duration && return (t=orientation*Inf,
            phi=orientation*a*Inf,
            tau=orientation*trapped_radial_tau)
        increment=_coordinate_increment(model,a,energy,lz,q,
            radial_value.radius,turn)
        return (t=orientation*increment.t,
            phi=orientation*increment.phi,
            tau=orientation*increment.tau)
    end
    domain=(mino=(-duration,duration),endpoint_closed=(true,true),
        endpoint_roles=(:past_horizon,:future_horizon),bl_mino=(-duration,duration))
    reference=(lambda0_event=:finite_turning_point,t_phi_zero_event=:finite_turning_point,
        t_phi_zero_lambda=0.0,t_phi_zero_radius=turn,tau_zero_event=:finite_turning_point,
        lambda_regular=duration)
    return (radial_r=radial_r, radial_sign=radial_sign, radial_history=radial_history,
            domain=domain, reference=reference,
            duration=duration, trapped_past_polar=trapped_past_polar,
            trapped_future_polar=trapped_future_polar, trapped_radial_tau=trapped_radial_tau)
end

# Horizon-endpoint record: regular-series matching for strict horizons, or the leading
# Laurent/Puiseux behaviour for repeated (double or triple) horizon roots.
function _extremal_endpoint_metadata(kind,endpoint,model,a,energy)
    return if endpoint!==nothing
        base=(kind=:exact_extremal_strict,
            series_order=endpoint.order,generated_order=endpoint.generated_order,
            v_match=endpoint.v_match,psi_match=endpoint.psi_match,
            matching_rule=:last_retained_and_orders_15_16_below_targets,
            regular_v_from_horizon=endpoint.v,
            regular_psi_from_horizon=endpoint.psi)
        kind===:trapped ? merge(base,(
            horizon_chart=(past=:retarded_u_chi,
                future=:advanced_v_psi),
            endpoint_zero=(past=(u=0.0,chi=0.0),
                future=(v=0.0,psi=0.0)),
            bl_endpoint_behavior=(:t_diverges_at_both_horizons,
                :phi_diverges_at_both_horizons))) : base
    elseif model!==nothing && model.horizon_root &&
            kind in (:horizon_root_turn,:horizon_root_from_infinity)
        c=model.h.h2; b=model.h.h1
        if abs(c)>1.0e-12
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
        polar_sector=nothing,polar_phase=0.0,axis=nothing,reference_radius=nothing)
    kind=spec.kind
    model=kind===:constant ? nothing :
        (kerr_geo_tier(spec.id) === :extremal ? _horizon_root_model(energy,q;
            turns=Tuple(x for x in (spec.lower,spec.upper) if isfinite(x) && x!=1)) :
         _strict_model(spec.id,a,energy,lz,q,structure;
            disposition=get(spec,:disposition,nothing)))
    polar=_extremal_polar_engine(a,energy,lz,q,polar_sector,polar_phase;axis=axis)
    has_strict_horizon=kind in (:horizon_to_turn,:horizon_to_repeated,
        :direct_capture,:trapped)
    endpoint=model===nothing || model.horizon_root || !has_strict_horizon ? nothing :
        _strict_endpoint_engine(model,a,energy,lz,q)


    builder = if kind===:constant
        _extremal_constant_branch
    elseif kind in (:periodic,:horizon_root_island)
        _extremal_periodic_branch
    elseif kind in (:scatter,:horizon_root_scatter,:homoclinic_outer,:horizon_root_turn)
        _extremal_two_sided_branch
    elseif kind in (:infinity_to_repeated,:horizon_root_from_infinity)
        _extremal_infall_branch
    elseif kind in (:horizon_to_turn,:horizon_to_repeated,:direct_capture)
        _extremal_horizon_branch
    elseif kind===:trapped
        _extremal_trapped_branch
    else
        error("Unknown exact-extremal motion kind $(kind).")
    end
    branch=builder(kind,spec,model,a,energy,lz,q,polar,reference_radius)
    (; radial_r,radial_sign,radial_history,domain)=branch
    reference=merge(branch.reference,(polar_phase=polar.metadata.phase,
        polar_phase_convention=polar.metadata.phase_convention))
    duration=get(branch,:duration,nothing)
    trapped_past_polar=get(branch,:trapped_past_polar,nothing)
    trapped_future_polar=get(branch,:trapped_future_polar,nothing)
    trapped_radial_tau=get(branch,:trapped_radial_tau,nothing)

    state_raw = function (lambda)
        local v, psi, rstar, tau      # not the trajectory closures of the same names below
        lam=_check_domain(lambda,domain)
        radius=radial_r(lam)
        pol=polar.formula(lam)
        if kind in (:horizon_to_turn,:horizon_to_repeated,:direct_capture)
            if lam==0
                return (t=Inf,r=1.0,z=pol.z,theta=pol.theta,
                    phi=a>0 ? Inf : -Inf,tau=0.0,rstar=-Inf,
                    v=0.0,psi=0.0,u=Inf,chi=a>0 ? Inf : -Inf)
            end
            v=-endpoint.v(radius)+pol.t
            psi=-endpoint.psi(radius)+pol.phi
            rstar=_rstar(radius); phi_h=_phi_h(a,radius)
            tau=-(model.basis(radius).I2-model.basis(1.0).I2)+pol.tau
            return (t=v-rstar,r=radius,z=pol.z,theta=pol.theta,
                phi=psi-phi_h,tau=tau,rstar=rstar,v=v,psi=psi,
                u=v-2rstar,chi=psi-2phi_h)
        end
        if kind===:trapped && abs(lam)==duration
            orientation=Base.sign(lam)
            tau_value=orientation*trapped_radial_tau+pol.tau
            if lam<0
                return (t=-Inf,r=1.0,z=pol.z,theta=pol.theta,
                    phi=-a*Inf,tau=tau_value,rstar=-Inf,
                    v=-Inf,psi=-a*Inf,u=0.0,chi=0.0)
            end
            return (t=Inf,r=1.0,z=pol.z,theta=pol.theta,
                phi=a*Inf,tau=tau_value,rstar=-Inf,
                v=0.0,psi=0.0,u=Inf,chi=a*Inf)
        end
        radial_value=radial_history(lam)
        time_value=radial_value.t+pol.t
        phi_value=radial_value.phi+pol.phi
        rstar_value=_rstar(radius); phi_h_value=_phi_h(a,radius)
        if kind===:trapped && lam<0
            u_value=endpoint.v(radius)+pol.t-trapped_past_polar.t
            chi_value=endpoint.psi(radius)+pol.phi-trapped_past_polar.phi
            return (t=time_value,r=radius,z=pol.z,theta=pol.theta,
                phi=phi_value,tau=radial_value.tau+pol.tau,
                rstar=rstar_value,v=time_value+rstar_value,
                psi=phi_value+phi_h_value,u=u_value,chi=chi_value)
        elseif kind===:trapped && lam>0
            v_value=-endpoint.v(radius)+pol.t-trapped_future_polar.t
            psi_value=-endpoint.psi(radius)+pol.phi-trapped_future_polar.phi
            return (t=time_value,r=radius,z=pol.z,theta=pol.theta,
                phi=phi_value,tau=radial_value.tau+pol.tau,
                rstar=rstar_value,v=v_value,psi=psi_value,
                u=time_value-rstar_value,chi=phi_value-phi_h_value)
        end
        return (t=time_value,r=radius,z=pol.z,theta=pol.theta,phi=phi_value,
            tau=radial_value.tau+pol.tau,rstar=rstar_value,
            v=time_value+rstar_value,psi=phi_value+phi_h_value,
            u=time_value-rstar_value,chi=phi_value-phi_h_value)
    end
    # Trapped orbits: v, psi are anchored on the future horizon and u, chi on the past
    # one. Each is a single function on the whole orbit (v - t - r_*, u - t + r_*, ... are
    # constants of the motion), so carry each anchor across the turning point.
    state = if kind===:trapped
        future=state_raw(duration/2); past=state_raw(-duration/2)
        cv=future.v-future.t-future.rstar; cpsi=future.psi-future.phi-_phi_h(a,future.r)
        cu=past.u-past.t+past.rstar; cpo=past.chi-past.phi+_phi_h(a,past.r)
        function (lambda)
            value=state_raw(lambda)
            lam=float(lambda)
            abs(lam)==duration && return value
            phi_h=_phi_h(a,value.r)
            lam<=0 && (value=merge(value,(v=value.t+value.rstar+cv,
                psi=value.phi+phi_h+cpsi)))
            lam>=0 && (value=merge(value,(u=value.t-value.rstar+cu,
                chi=value.phi-phi_h+cpo)))
            return value
        end
    else
        state_raw
    end

    r=lambda -> radial_r(_check_domain(lambda,domain))
    z=lambda -> polar.formula(_check_domain(lambda,domain)).z
    theta=lambda -> polar.formula(_check_domain(lambda,domain)).theta
    t=lambda -> state(lambda).t
    phi=lambda -> state(lambda).phi
    tau=lambda -> state(lambda).tau
    rstar=lambda -> state(lambda).rstar
    v=lambda -> state(lambda).v
    psi=lambda -> state(lambda).psi
    u=lambda -> state(lambda).u
    chi=lambda -> state(lambda).chi
    kin=_kinematics(a,energy,lz,q;r=r,z=z,
        uz=lambda -> polar.formula(_check_domain(lambda,domain)).uz,
        sin2=lambda -> polar.formula(_check_domain(lambda,domain)).sin2,
        sign_r=lambda -> radial_sign(_check_domain(lambda,domain)),
        R=radius -> kerr_radial_potential(a,energy,lz,q,radius))

    endpoint_metadata=_extremal_endpoint_metadata(kind,endpoint,model,a,energy)
    roots_metadata=(radial=structure===nothing ? () : structure.real_roots,
        raw=structure===nothing ? () : structure.raw_roots,
        lower=spec.lower,upper=spec.upper,polar=(sector=polar.sector,))
    metric_limit=a==1 ? :extremal_plus : :extremal_minus
    (;velocity,potentials,residuals)=_kinematic_fields(kin)
    return _member(spec.broad,spec.id;tier=:extremal,
        constants=(a=a,E=energy,Lz=lz,Q=q),roots=roots_metadata,reference=reference,
        domain=domain,
        trajectory=(t=t,r=r,theta=theta,z=z,phi=phi,tau=tau,rstar=rstar,v=v,psi=psi,u=u,
            chi=chi),
        velocity,potentials,residuals,
        status=(supported=true,metric_limit=metric_limit,
            metric_dispatch=a==1 ? :exact_positive : :exact_negative,
            formula_family=spec.formula,motion_kind=kind,endpoint=endpoint_metadata,
            disposition_id=get(spec,:disposition,nothing),
            polar_sector=polar.sector),
        spectral=SpectralStatus(() -> (_polar_spectral(polar),)))
end

"""
    kerr_geo_extremal_family(a, E, Lz, Q; polar_sector=nothing, polar_phase=0.0,
                             axis=nothing, reference_radius=nothing)

Classify `(E, Lz, Q)` at exact `a = +1` or `a = -1` and construct every admitted member.
With `P_H = 2E − aLz = 0` the members are the horizon-root cases A-X1, A-X2, B-X1, B-X2,
C-X1…C-X4, D-X1, D-X2; otherwise they keep their primary case IDs. `polar_sector` selects the
polar sector, `polar_phase` is the polar phase at `λ = 0`, `axis = :north` or `:south` puts
the motion on the spin axis (`Lz = 0`, `Q = a²(1 − E²)`), and `reference_radius` is the
radius at `λ = 0` for members that come in from infinity to a repeated or horizon root
(K7, K10, C-X1…C-X4).
"""
function kerr_geo_extremal_family(a::Real,energy::Real,lz::Real,q::Real;
        polar_sector=nothing,polar_phase=0.0,axis=nothing,reference_radius=nothing)
    (a==1 || a==-1) || throw(DomainError(a,
        "Exact-extremal family dispatch requires exact a=+1 or a=-1."))
    if axis!==nothing
        abs(lz)<=1e-12 || error("Axis initial data requires Lz=0.")
        qaxis=kerr_axis_carter_q(a,energy)
        abs(q-qaxis)<=1e-10*max(1.0,abs(qaxis)) || error(
            "Axis initial data requires Q=a^2(1-E^2).")
    end
    classification=_classification(a,float(energy),float(lz),float(q))
    members=Tuple(_make_trajectory(float(a),float(energy),float(lz),float(q),
        spec,classification.structure;polar_sector=polar_sector,
        polar_phase=polar_phase,axis=axis,reference_radius=reference_radius)
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
