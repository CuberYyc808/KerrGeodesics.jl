# Polar motion of E ≥ 1 scatter and capture orbits: `_window_polar_motion` (the pendular
# polar-engine solution with the E = 1 phase convention) and the hyperbolic (E > 1) and
# parabolic (E = 1 or a = 0) closed forms used by the finite-window scatter construction.


# The polar phase at the reference radial event and the polar event at phase zero, taken
# from `_window_polar_motion` (northern turning point for E > 1, northward equator crossing
# at E = 1, `:not_applicable` for equatorial orbits).
_polar_phase_metadata(polar, polar_phase, reference_event) =
    (polar_phase=polar.inclined ? polar_phase : nothing, reference_radial_event=reference_event,
     phase_convention=polar.metadata.phase_convention)

# E > 1 polar roots: z_+ in (0,1) (turning point) and z_- < 0, with
# m = z_+/(z_+ - z_-) and omega = sqrt(beta (z_+ - z_-)), beta = a^2(E^2 - 1) = -c;
# formed from c*u so that nothing is divided by the small beta.
function _hyperbolic_polar_roots(a, energy, lz, q)
    roots = _polar_quadratic_roots(a, energy, lz, q)
    roots.disc >= 0 || error("Generic inclined hyperbolic polar roots are not real.")
    pairs = ((roots.u_small, roots.cu_small), (roots.u_big, roots.cu_big))
    (zplus, cplus), (zminus, cminus) = pairs[1][1] >= pairs[2][1] ? pairs : reverse(pairs)
    # (z_+ = 1 on the axis, Lz = 0, or when Lz² ≪ Q rounds it there; the polar engine
    # forms 1 − z_+ separately with _polar_one_minus_root)
    0 < zplus <= 1 || error("Generic inclined hyperbolic polar root z_+ is outside (0,1].")
    zminus < 0 || error("Generic inclined hyperbolic polar root z_- must be negative.")
    spread = cminus - cplus                  # beta (z_+ - z_-) = -c z_+ + c z_-
    return (zplus=zplus, zminus=zminus, m=-cplus / spread, omega=sqrt(spread))
end

function _hyperbolic_polar_parameters(a, energy, lz, q)
    beta = a^2 * _e2m1(energy)
    beta > 0 || error("Generic inclined hyperbolic polar motion requires a^2(E^2-1)>0.")
    (; zplus, zminus, m, omega) = _hyperbolic_polar_roots(a, energy, lz, q)
    kcomplete = Elliptic.F(pi / 2, m)
    ecomplete = Elliptic.E(pi / 2, m)
    nphi = -zplus / (1 - zplus)
    picomplete = Elliptic.Pi(nphi, pi / 2, m)
    period_u = 2 * kcomplete
    phi_period = 2 * lz * picomplete / (omega * (1 - zplus))
    t_period = 2 / omega * ((a * lz - a^2 * energy) * kcomplete +
                             a^2 * energy * zplus / m *
                             (ecomplete + (m - 1) * kcomplete))
    tau_period = 2 * a^2 * zplus / (omega * m) *
        (ecomplete + (m - 1) * kcomplete)
    return (
        beta=beta,
        zplus=zplus,
        zminus=zminus,
        m=m,
        omega=omega,
        kcomplete=kcomplete,
        ecomplete=ecomplete,
        nphi=nphi,
        picomplete=picomplete,
        period_u=period_u,
        phi_period=phi_period,
        t_period=t_period,
        tau_period=tau_period,
    )
end

function _hyperbolic_amplitude_from_quarter_u(u, m; max_iter=70)
    lo = 0.0
    hi = pi / 2
    for _ in 1:max_iter
        mid = 0.5 * (lo + hi)
        if Elliptic.F(mid, m) < u
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

function _hyperbolic_jacobi_sncndn_from_u(u, m, kcomplete)
    period = 4 * kcomplete
    x = mod(u, period)
    if x <= kcomplete
        phi = _hyperbolic_amplitude_from_quarter_u(x, m)
        sn = sin(phi)
        cn = cos(phi)
    elseif x <= 2 * kcomplete
        phi = _hyperbolic_amplitude_from_quarter_u(2 * kcomplete - x, m)
        sn = sin(phi)
        cn = -cos(phi)
    elseif x <= 3 * kcomplete
        phi = _hyperbolic_amplitude_from_quarter_u(x - 2 * kcomplete, m)
        sn = -sin(phi)
        cn = -cos(phi)
    else
        phi = _hyperbolic_amplitude_from_quarter_u(4 * kcomplete - x, m)
        sn = -sin(phi)
        cn = cos(phi)
    end
    dn = sqrt(max(1 - m * sn^2, 0.0))
    return sn, cn, dn
end

function _hyperbolic_polar_low_half_primitive(a, energy, lz, pars, y)
    yy = clamp(y, 0.0, pars.kcomplete)
    psi = _hyperbolic_amplitude_from_quarter_u(yy, pars.m)
    f = Elliptic.F(psi, pars.m)
    e = Elliptic.E(psi, pars.m)
    piinc = Elliptic.Pi(pars.nphi, psi, pars.m)
    phi = lz * piinc / (pars.omega * (1 - pars.zplus))
    t = ((a * lz - a^2 * energy) * f +
         a^2 * energy * pars.zplus / pars.m *
         (e + (pars.m - 1) * f)) / pars.omega
    tau = a^2 * pars.zplus / pars.m *
        (e + (pars.m - 1) * f) / pars.omega
    return (phi=phi, t=t, tau=tau)
end

function _hyperbolic_polar_reduced_primitive(a, energy, lz, pars, y)
    yy = y
    if yy < 0.0 && yy > -2e-14
        yy = 0.0
    elseif yy > pars.period_u && yy < pars.period_u + 2e-14
        yy = pars.period_u
    end
    0.0 <= yy <= pars.period_u + 1e-12 ||
        error("Polar reduced argument lies outside [0,2K].")
    if yy <= pars.kcomplete
        return _hyperbolic_polar_low_half_primitive(a, energy, lz, pars, yy)
    end
    reflected = _hyperbolic_polar_low_half_primitive(
        a,
        energy,
        lz,
        pars,
        pars.period_u - yy,
    )
    return (phi=pars.phi_period - reflected.phi,
        t=pars.t_period - reflected.t,
        tau=pars.tau_period - reflected.tau)
end

function _hyperbolic_polar_global_primitive(a, energy, lz, pars, u)
    cycle = floor(Int, u / pars.period_u)
    y = u - cycle * pars.period_u
    if y < 0.0
        cycle -= 1
        y += pars.period_u
    end
    reduced = _hyperbolic_polar_reduced_primitive(a, energy, lz, pars, y)
    return (
        phi=cycle * pars.phi_period + reduced.phi,
        t=cycle * pars.t_period + reduced.t,
        tau=cycle * pars.tau_period + reduced.tau,
    )
end

function _hyperbolic_polar_formula(a, energy, lz, pars, u0, start_primitive, lambda)
    u = u0 + pars.omega * lambda
    sn, cn, dn = _hyperbolic_jacobi_sncndn_from_u(u, pars.m, pars.kcomplete)
    z = sqrt(pars.zplus) * cn
    uz = -sqrt(pars.zplus) * pars.omega * sn * dn
    primitive = _hyperbolic_polar_global_primitive(a, energy, lz, pars, u)
    return (
        z=z,
        uz=uz,
        theta=acos(clamp(z, -1.0, 1.0)),
        phi=primitive.phi - start_primitive.phi,
        t=primitive.t - start_primitive.t,
        tau=primitive.tau - start_primitive.tau,
    )
end

function _parabolic_polar_parameters(lz, q)
    q > 0 || error("Generic inclined parabolic polar motion requires Q>0.")
    omega = sqrt(q + lz^2)
    return (omega=omega, qamp=sqrt(q / (q + lz^2)), c=abs(lz) / omega)
end

function _parabolic_phi_angle_primitive(lz, pars, u)
    abs(lz) <= 1e-14 && return 0.0
    cycle = floor(Int, u / pi)
    y = u - cycle * pi
    if y < 0.0
        cycle -= 1
        y += pi
    end
    local_angle = if y <= pi / 2
        atan(pars.c * tan(y))
    else
        pi - atan(pars.c * tan(pi - y))
    end
    return sign(lz) * (cycle * pi + local_angle)
end

function _parabolic_t_angle_primitive(a, lz, pars, u)
    qamp2 = pars.qamp^2
    return ((a * lz - a^2) * u +
            0.5 * a^2 * qamp2 * (u - sin(u) * cos(u))) / pars.omega
end

function _parabolic_tau_angle_primitive(a, pars, u)
    return 0.5 * a^2 * pars.qamp^2 *
        (u - sin(u) * cos(u)) / pars.omega
end

function _parabolic_polar_formula(a, lz, pars, u0, start_phi, start_t, lambda)
    u = u0 + pars.omega * lambda
    z = pars.qamp * sin(u)
    uz = pars.qamp * pars.omega * cos(u)
    return (
        z=z,
        uz=uz,
        theta=acos(clamp(z, -1.0, 1.0)),
        phi=_parabolic_phi_angle_primitive(lz, pars, u) - start_phi,
        t=_parabolic_t_angle_primitive(a, lz, pars, u) - start_t,
        tau=_parabolic_tau_angle_primitive(a, pars, u) -
            _parabolic_tau_angle_primitive(a, pars, u0),
    )
end

"""
    _window_polar_motion(a, energy, lz, q, polar_phase)

Polar part of an E ≥ 1 orbit as Mino-time closures. Equatorial orbits (`Q = 0`) stay at
θ = π/2. Inclined orbits use the pendular polar-engine solution with phase 0 at the northern
turning point; at E = 1 the phase is shifted so that z = √z₋ sin(phase + ωλ), phase 0 on the
equator moving north. `polar_phase` is the phase at λ = 0.
`rootdata` is the polar entry of the orbit's `Roots` record.
"""
function _window_polar_motion(a, energy, lz, q, polar_phase)
    if abs(q) <= 1e-12
        return (inclined=false, z=λ -> 0.0, uz=λ -> 0.0, theta=λ -> pi / 2,
                phi=λ -> lz * λ, t=λ -> (a * lz - a^2 * energy) * λ, tau=λ -> 0.0,
                phidot=λ -> lz, tdot=λ -> a * lz - a^2 * energy, rootdata=nothing,
                primitive=λ -> ((a * lz - a^2 * energy) * λ, lz * λ, 0.0),
                position=λ -> (0.0, 0.0, 1.0),
                metadata=(sector=:equatorial, phase=0.0, phase_convention=:not_applicable))
    end
    # polar engine (PolarEngine.jl): z = √z₊ cn(phase + ωλ) for E > 1, phase 0 at the
    # northern turning point; at E = 1 the phase is shifted by K so that
    # z = √z₋ sin(phase + ωλ), phase 0 on the equator moving north
    sol = _polar_solution(a, energy, lz, q, :pendular, float(polar_phase))
    if abs(energy - 1) <= 1e-12
        sol = _polar_solution(a, energy, lz, q, :pendular,
            float(polar_phase) - Elliptic.K(sol.metadata.modulus))
    end
    md = sol.metadata
    beta_zero = abs(a^2 * _e2m1(energy)) <= 1e-14
    rootdata = (zplus=beta_zero ? md.zminus : md.zplus, zminus=beta_zero ? 0.0 : md.zminus,
        modulus=md.modulus, omega=md.omega, beta_zero=beta_zero)
    position = sol.position
    primitive = sol.primitive
    zf(λ) = position(λ)[1]
    return (inclined=true, z=zf, uz=λ -> position(λ)[2],
            theta=λ -> acos(clamp(zf(λ), -1.0, 1.0)), phi=λ -> primitive(λ)[2],
            t=λ -> primitive(λ)[1], tau=λ -> primitive(λ)[3],
            phidot=λ -> lz / position(λ)[3],
            tdot=λ -> a * lz - a^2 * energy * position(λ)[3],
            rootdata=rootdata, primitive=primitive, position=position,
            metadata=merge(md, (phase=float(polar_phase),
                phase_convention=abs(energy - 1) <= 1e-12 ?
                    :equator_northward_at_zero_phase : :northern_turning_point_at_zero_phase)))
end
