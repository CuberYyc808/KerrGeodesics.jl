# Polar motion of every sector (pendular, vortical, equatorial, axis crossing, constant
# latitude, equator attractive) as polar-engine solutions, and the choice of sector.

function _constant_latitude_polar_solution(a, energy, lz, q, polar_phase;
        hemisphere=:north)
    hemisphere in (:north, :south) || error(
        "Constant-latitude motion requires hemisphere=:north or :south.")
    beta = a^2 * _e2m1(energy)
    beta > 0 || error("Non-equatorial constant latitude requires E>1 and a!=0.")
    geometry = _polar_vortical_geometry(a, energy, lz, q)
    geometry.repeated || error(
        "Constants do not satisfy the constant-latitude double-root condition.")
    u0 = 0.5 * geometry.root_sum
    0 < u0 < 1 || error("Constant latitude must satisfy 0<u0<1.")
    signz = hemisphere === :north ? 1.0 : -1.0
    z0 = signz * sqrt(u0)
    phi_rate = lz / (1 - u0)
    time_rate = a * lz - a^2 * energy * (1 - u0)
    formula = lambda -> (
        z=z0,
        uz=0.0,
        sin2=1 - u0,
        theta=atan(sqrt(1 - u0), z0),
        phi=phi_rate * float(lambda),
        t=time_rate * float(lambda),
        tau=a^2 * u0 * float(lambda),
    )
    return (formula=formula, metadata=(
        sector=:constant_latitude,
        reading=geometry.reading,
        hemisphere=hemisphere,
        u0=u0,
        phase=polar_phase,
        phase_convention=:not_applicable,
        formula_kind=:double_polar_root,
    ))
end

function _vortical_polar_solution(a, energy, lz, q, polar_phase;
        hemisphere=:north)
    hemisphere in (:north, :south) || error(
        "Vortical motion requires hemisphere=:north or :south.")
    beta = a^2 * _e2m1(energy)
    beta > 0 && q < 0 || error("Vortical motion requires E>1 and Q<0.")
    roots = _polar_quadratic_roots(a, energy, lz, q)
    roots.disc > 0 || return _constant_latitude_polar_solution(
        a, energy, lz, q, polar_phase; hemisphere=hemisphere)
    uminus, uplus = minmax(roots.u_small, roots.u_big)
    # Θ(1) = −Lz² puts u₊ ≤ 1; for Lz² ≪ |Q| the rounded root can land an ulp above 1
    0 < uminus < uplus <= 1 + POLAR_ROOT_SLACK || error(
        "Vortical polar roots must satisfy 0<uminus<uplus≤1.")
    uplus = min(uplus, 1.0)
    # z = ±√u₊ dn(u|m), m = (u₊ − u₋)/u₊ with u₊ − u₋ = √disc/β, and 1 − m = u₋/u₊
    return _elliptic_polar_solution(a, energy, lz, q; kind=:dn, A=uplus,
        one_minus_A=_polar_one_minus_root(a, energy, lz, q, uplus),
        m=sqrt(roots.disc) / (beta * uplus), m1=uminus / uplus, omega=sqrt(beta * uplus),
        u0=float(polar_phase), sign=hemisphere === :north ? 1.0 : -1.0,
        metadata=(sector=:vortical, hemisphere=hemisphere, uminus=uminus,
            uplus=uplus, phase=polar_phase,
            phase_convention=:maximum_absolute_latitude_at_zero_phase))
end

function _equator_attractive_polar_solution(a, energy, lz, q, polar_phase;
        hemisphere=:north)
    hemisphere in (:north, :south) || error(
        "Equator-attractive motion requires hemisphere=:north or :south.")
    beta = a^2 * _e2m1(energy)
    beta > 0 && iszero(q) || error(
        "Equator-attractive motion requires E>1 and Q=0.")
    !iszero(lz) || error(
        "The Lz=0 limit crosses the axis and uses the axis-crossing formula.")
    lz^2 < beta || error("Equator-attractive constants require 0<Lz^2<a^2(E^2-1).")
    # the minimum of sin²θ, ε = Lz²/β, is kept from the constants: 1 − (1 − ε) would lose it
    complement = lz^2 / beta
    amplitude2 = 1 - complement
    omega = sqrt(beta * amplitude2)
    ratio = sqrt(amplitude2 * beta) / abs(lz)             # √((1 − ε)/ε)
    signz = hemisphere === :north ? 1.0 : -1.0
    function primitive(u)
        y = tanh(u)
        phi = lz / omega * u + sign(lz) * atan(y * ratio)
        time = ((a * lz - a^2 * energy) * u +
            a^2 * energy * amplitude2 * y) / omega
        tau = a^2 * amplitude2 * y / omega
        return (phi=phi, t=time, tau=tau)
    end
    u0 = float(polar_phase)
    p0 = primitive(u0)
    formula = function(lambda)
        u = u0 + omega * float(lambda)
        sech = inv(cosh(u))
        z = signz * sqrt(amplitude2) * sech
        uz = -signz * sqrt(amplitude2) * omega * sech * tanh(u)
        p = primitive(u)
        sin2=complement + amplitude2 * tanh(u)^2
        return (z=z, uz=uz, sin2=sin2,
            theta=atan(sqrt(sin2), z), phi=p.phi - p0.phi, t=p.t - p0.t,
            tau=p.tau - p0.tau)
    end
    return (formula=formula, metadata=(
        sector=:equator_attractive,
        hemisphere=hemisphere,
        amplitude_squared=amplitude2,
        omega=omega,
        phase=polar_phase,
        phase_convention=:maximum_absolute_latitude_at_zero_phase,
        formula_kind=:hyperbolic_sech,
    ))
end

function _axis_crossing_polar_solution(a, energy, lz, q, polar_phase;
        hemisphere=:north)
    _zero_lz(a, energy, lz, q) || error(
        "Axis-crossing polar motion requires Lz=0.")
    hemisphere in (:north, :south) || error(
        "Axis-crossing polar motion requires hemisphere=:north or :south.")
    meta = (sector=:axis_crossing, hemisphere=hemisphere, phase=polar_phase,
        phase_convention=:north_axis_at_zero_phase,
        coordinate_chart=:z_primary_theta_folded_phi_patch_required_at_axis)
    if !iszero(lz)
        # |Lz| inside the axis tolerance but not zero: the orbit turns just short of the
        # pole. The generic forms keep 1 − z_turn exactly, hence the Lz/(1 − z²) spike.
        sol = q < 0 ?
            _vortical_polar_solution(a, energy, lz, q, polar_phase; hemisphere=hemisphere) :
            q > 0 ? _polar_solution(a, energy, lz, q, :pendular, polar_phase) :
            _equator_attractive_polar_solution(a, energy, lz, q, polar_phase;
                hemisphere=hemisphere)
        return merge(sol, (metadata=merge(sol.metadata, meta),))
    end
    c = -a^2 * _e2m1(energy)                      # Θ = q(1 − z²) − c z²(1 − z²) for Lz = 0
    if c >= 0
        q >= c || error("Lz=0 axis crossing with E<=1 requires Q>=a^2(1-E^2).")
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:cd, A=1.0,
            one_minus_A=0.0, m=c / q, m1=(q - c) / q, omega=sqrt(q), u0=float(polar_phase),
            metadata=meta)
    end
    beta = -c
    beta + q > 0 || error("Axis-crossing motion with E>1 requires Q>-a^2(E^2-1).")
    if q > 0
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:cn, A=1.0,
            one_minus_A=0.0, m=beta / (beta + q), m1=q / (beta + q), omega=sqrt(beta + q),
            complementary_modulus=sqrt(q)/sqrt(beta+q),
            u0=float(polar_phase), metadata=meta)
    elseif q < 0
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:dn, A=1.0,
            one_minus_A=0.0, m=(beta + q) / beta, m1=-q / beta, omega=sqrt(beta),
            u0=float(polar_phase), sign=hemisphere === :north ? 1.0 : -1.0,
            metadata=meta)
    end
    # q = 0: z = ±sech(u), elementary primitives
    omega = sqrt(beta)
    u0 = float(polar_phase)
    signz = hemisphere === :north ? 1.0 : -1.0
    primitive(u) = (t=-a^2 * energy * (u - tanh(u)) / omega, tau=a^2 * tanh(u) / omega)
    p0 = primitive(u0)
    formula = function(lambda)
        u = u0 + omega * float(lambda)
        sech = inv(cosh(u))
        z = signz * sech
        p = primitive(u)
        sin2=tanh(u)^2
        return (z=z, uz=-signz * omega * sech * tanh(u), sin2=sin2,
            theta=atan(sqrt(sin2), z), phi=0.0,
            t=p.t - p0.t, tau=p.tau - p0.tau)
    end
    return (formula=formula, metadata=merge(meta, (modulus=1.0, omega=omega,
        formula_kind=:hyperbolic_sech_axis_crossing)))
end

# The polar sector of a member: the requested one, else the classified one; when the classifier
# needs initial data to choose, the sector the constants allow on their own.
function _resolved_polar_sector(classification, requested, energy, lz, q)
    requested !== nothing && return requested
    selected = classification.PolarMetadata.selected
    return selected !== :polar_initial_data_required ? selected :
        _constants_polar_sector(energy, lz, q)
end

# equatorial for Q = 0 (Lz ≠ 0), pendular for Q > 0, vortical for Q < 0 and E > 1: the sign
# of Q decides, and only Q = 0 itself is the equatorial limit
_constants_polar_sector(energy, lz, q) =
    iszero(q) && !iszero(lz) ? :equatorial : q > 0 ? :pendular :
    q < 0 && energy > 1 ? :vortical : nothing

function _polar_solution(a, energy, lz, q, sector, polar_phase;
        hemisphere=:north)
    sector === :polar_initial_data_required && error(
        "These constants admit several polar sectors; choose one with the polar_sector keyword.")
    sector === :axis_constant && error(
        "Axis trajectories use the separate axis-infall constructor.")
    sector === :axis_crossing && return _axis_crossing_polar_solution(
        a, energy, lz, q, polar_phase; hemisphere=hemisphere)
    sector === :equator_attractive && return _equator_attractive_polar_solution(
        a, energy, lz, q, polar_phase; hemisphere=hemisphere)
    sector === :constant_latitude && return _constant_latitude_polar_solution(
        a, energy, lz, q, polar_phase; hemisphere=hemisphere)
    if sector === :vortical
        beta = a^2 * _e2m1(energy)
        geometry = _polar_vortical_geometry(a, energy, lz, q)
        return geometry.repeated ?
            _constant_latitude_polar_solution(
                a, energy, lz, q, polar_phase; hemisphere=hemisphere) :
            _vortical_polar_solution(
                a, energy, lz, q, polar_phase; hemisphere=hemisphere)
    end
    sector === :equatorial && return _equatorial_polar_solution(a, energy, lz)
    q > 0 || error("Pendular polar motion requires Q>0.")

    c = -a^2 * _e2m1(energy)
    if c >= 0
        # z = √z₋ cd(u|m): E < 1 (c > 0) and parabolic / Schwarzschild (c = 0, m = 0)
        roots = _polar_quadratic_roots(a, energy, lz, q)
        roots.disc >= 0 || error("Polar roots with E<=1 are not real.")
        # (z₋ = 1 up to rounding when Lz² ≪ Q; 1 − z₋ is formed separately)
        zminus = roots.u_small
        # The amplitude remains representable when its square underflows.
        amplitude = sqrt(q) / sqrt(roots.cu_big)
        0 < amplitude <= sqrt(1 + POLAR_ROOT_SLACK) || error("Polar turning root with E<=1 lies outside (0,1].")
        zminus = min(zminus, 1.0)
        # m = z₋/z₊ and 1 − m = (z₊ − z₋)/z₊ = √disc/(c z₊)
        return _elliptic_polar_solution(a, energy, lz, q; kind=:cd, A=zminus,
            amplitude=min(amplitude, 1.0),
            one_minus_A=_polar_one_minus_root(a, energy, lz, q, zminus),
            m=roots.cu_small / roots.cu_big, m1=sqrt(roots.disc) / roots.cu_big,
            omega=sqrt(roots.cu_big), u0=float(polar_phase),
            metadata=(sector=:pendular, zminus=zminus, zplus=roots.u_big, phase=polar_phase,
                phase_convention=:northern_turning_point_at_zero_phase))
    end
    # z = √z₊ cn(u|m): E > 1
    (; zplus, zminus, m, m1, omega, amplitude, complementary_modulus) =
        _hyperbolic_polar_roots(a, energy, lz, q)
    return _elliptic_polar_solution(a, energy, lz, q; kind=:cn, A=zplus,
        amplitude, complementary_modulus,
        one_minus_A=_polar_one_minus_root(a, energy, lz, q, zplus), m=m, m1=m1,
        omega=omega, u0=float(polar_phase),
        metadata=(sector=:pendular, zminus=zminus, zplus=zplus,
            phase=polar_phase, phase_convention=:northern_turning_point_at_zero_phase))
end
