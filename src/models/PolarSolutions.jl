# Polar motion of every sector (pendular, vortical, equatorial, axis crossing, constant
# latitude, equator attractive) as polar-engine solutions, and the choice of sector.

function _constant_latitude_polar_solution(a, energy, lz, q, polar_phase;
        hemisphere=:north)
    hemisphere in (:north, :south) || error(
        "Constant-latitude motion requires hemisphere=:north or :south.")
    beta = a^2 * _e2m1(energy)
    beta > 0 || error("Non-equatorial constant latitude requires E>1 and a!=0.")
    root_sum = (beta - q - lz^2) / beta
    root_product = -q / beta
    discriminant = root_sum^2 - 4 * root_product
    tolerance = 1.0e-10 * max(1.0, root_sum^2, abs(root_product))
    abs(discriminant) <= tolerance || error(
        "Constants do not satisfy the constant-latitude double-root condition.")
    u0 = 0.5 * root_sum
    0 < u0 < 1 || error("Constant latitude must satisfy 0<u0<1.")
    expected_lz2 = beta * (1 - u0)^2
    expected_q = -beta * u0^2
    abs(lz^2 - expected_lz2) <= 1.0e-9 * max(1.0, expected_lz2) ||
        error("Constant-latitude Lz condition failed.")
    abs(q - expected_q) <= 1.0e-9 * max(1.0, abs(expected_q)) ||
        error("Constant-latitude Q condition failed.")
    signz = hemisphere === :north ? 1.0 : -1.0
    z0 = signz * sqrt(u0)
    phi_rate = lz / (1 - u0)
    time_rate = a * lz - a^2 * energy * (1 - u0)
    formula = lambda -> (
        z=z0,
        uz=0.0,
        sin2=1 - u0,
        theta=acos(z0),
        phi=phi_rate * float(lambda),
        t=time_rate * float(lambda),
        tau=a^2 * u0 * float(lambda),
    )
    return (formula=formula, metadata=(
        sector=:constant_latitude,
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
    0 < uminus < uplus <= 1 || error(
        "Vortical polar roots must satisfy 0<uminus<uplus≤1.")
    # z = ±√u₊ dn(u|m), m = (u₊ − u₋)/u₊, with u₊ − u₋ = √disc/β
    return _elliptic_polar_solution(a, energy, lz, q; kind=:dn, A=uplus,
        one_minus_A=_polar_one_minus_root(a, energy, lz, q, uplus),
        m=sqrt(roots.disc) / (beta * uplus), omega=sqrt(beta * uplus),
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
    beta > 0 && abs(q) <= 1.0e-13 || error(
        "Equator-attractive motion requires E>1 and Q=0.")
    abs(lz) > 1.0e-13 || error(
        "The Lz=0 limit crosses the axis and uses the axis-crossing formula.")
    amplitude2 = 1 - lz^2 / beta
    0 < amplitude2 < 1 || error(
        "Equator-attractive constants require 0<1-Lz^2/beta<1.")
    complement = 1 - amplitude2
    omega = sqrt(beta * amplitude2)
    signz = hemisphere === :north ? 1.0 : -1.0
    function primitive(u)
        y = tanh(u)
        phi = lz / omega * (u + sqrt(amplitude2 / complement) *
            atan(y * sqrt(amplitude2 / complement)))
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
        return (z=z, uz=uz, sin2=complement + amplitude2 * tanh(u)^2,
            theta=acos(clamp(z, -1.0, 1.0)), phi=p.phi - p0.phi, t=p.t - p0.t,
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
    abs(lz) <= _axis_lz_tolerance(q, -a^2 * _e2m1(energy)) || error(
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
            _polar_solution(a, energy, lz, q, :pendular, polar_phase)
        return merge(sol, (metadata=merge(sol.metadata, meta),))
    end
    c = -a^2 * _e2m1(energy)                      # Θ = q(1 − z²) − c z²(1 − z²) for Lz = 0
    if c >= 0
        q >= c || error("Lz=0 axis crossing with E<=1 requires Q>=a^2(1-E^2).")
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:cd, A=1.0,
            one_minus_A=0.0, m=c / q, omega=sqrt(q), u0=float(polar_phase), metadata=meta)
    end
    beta = -c
    beta + q > 0 || error("Axis-crossing motion with E>1 requires Q>-a^2(E^2-1).")
    if q > 0
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:cn, A=1.0,
            one_minus_A=0.0, m=beta / (beta + q), omega=sqrt(beta + q),
            u0=float(polar_phase), metadata=meta)
    elseif q < 0
        return _elliptic_polar_solution(a, energy, 0.0, q; kind=:dn, A=1.0,
            one_minus_A=0.0, m=(beta + q) / beta, omega=sqrt(beta),
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
        return (z=z, uz=-signz * omega * sech * tanh(u), sin2=tanh(u)^2,
            theta=acos(clamp(z, -1.0, 1.0)), phi=0.0,
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

# equatorial for Q = 0 (Lz ≠ 0), pendular for Q > 0, vortical for Q < 0 and E > 1
_constants_polar_sector(energy, lz, q) =
    abs(q) <= 1.0e-12 && abs(lz) > 1.0e-12 ? :equatorial : q > 0 ? :pendular :
    q < 0 && energy > 1 ? :vortical : nothing

function _polar_solution(a, energy, lz, q, sector, polar_phase;
        hemisphere=:north)
    abs(q) <= 1.0e-13 && abs(lz) > 1.0e-13 && sector !== :equator_attractive &&
        return _equatorial_polar_solution(a, energy, lz)
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
        root_sum = (beta - q - lz^2) / beta
        root_product = -q / beta
        discriminant = root_sum^2 - 4 * root_product
        tolerance = 1.0e-10 * max(1.0, root_sum^2, abs(root_product))
        return abs(discriminant) <= tolerance ?
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
        0 < zminus <= 1 + 4eps() || error("Polar turning root with E<=1 lies outside (0,1].")
        zminus = min(zminus, 1.0)
        m = roots.cu_small / roots.cu_big
        return _elliptic_polar_solution(a, energy, lz, q; kind=:cd, A=zminus,
            one_minus_A=_polar_one_minus_root(a, energy, lz, q, zminus), m=m,
            omega=sqrt(roots.cu_big), u0=float(polar_phase),
            metadata=(sector=:pendular, zminus=zminus, zplus=roots.u_big, phase=polar_phase,
                phase_convention=:northern_turning_point_at_zero_phase))
    end
    # z = √z₊ cn(u|m): E > 1
    (; zplus, zminus, m, omega) = _hyperbolic_polar_roots(a, energy, lz, q)
    return _elliptic_polar_solution(a, energy, lz, q; kind=:cn, A=zplus,
        one_minus_A=_polar_one_minus_root(a, energy, lz, q, zplus), m=m,
        omega=omega, u0=float(polar_phase),
        metadata=(sector=:pendular, zminus=zminus, zplus=zplus,
            phase=polar_phase, phase_convention=:northern_turning_point_at_zero_phase))
end

