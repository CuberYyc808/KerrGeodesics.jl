# Finite-window scatter API `kerr_geo_scatter` (E >= 1, infinity -> turning point ->
# infinity), its closed-form radial pieces and the asymptotic diagnostics.

"""
    KerrGeoScatter

E ≥ 1 orbit coming in from infinity, turning at closest approach (λ = 0) and
escaping again. `Formula` is `:hyperbolic_scatter` (E > 1, four real roots) or
`:parabolic_scatter` (E = 1, three real roots). `kerr_geo_scatter_asymptotic_diagnostics`
gives the asymptotic directions and deflection angles.
"""
struct KerrGeoScatter
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

function _scatter_domain(lambda_infinity)
    return (mino=(-lambda_infinity, lambda_infinity),
            endpoint_roles=(:past_infinity, :future_infinity), endpoint_closed=(false, false))
end

function _d2_roots(outcome)
    roots = collect(outcome.roots)
    length(roots) == 4 || return nothing
    sort!(roots)
    return (rA=roots[1], rB=roots[2], rC=roots[3], rD=roots[4])
end

# Mino time from the turning point of a scattering leg, in [0, λ∞): within MINO_ENDPOINT_TOL
# below 0 it is the turning point
function _leg_time(s, λ∞)
    -MINO_ENDPOINT_TOL <= s < λ∞ || throw(DomainError(s,
        "Mino time from the turning point must lie in [0, λ∞ = $(λ∞))."))
    return max(s, 0.0)
end

function _d2_rminus(a)
    return 1 - sqrt(max(0.0, 1 - a^2))
end

# the pole primitive of 1/(r − h) on the D2 leg: F and the R_J term share the common R_F, no
# division by r_C − h remains, and cos φ is not recovered from a sine close to one
function _d2_pole_primitive(leg, h, phi)
    roots = leg.roots
    n1 = _four_real_characteristic(leg, h)[2]
    s, c = sincos(phi)
    alpha = (roots.rD - roots.rC) * (roots.rD - roots.rA) / ((roots.rC - roots.rA) * (roots.rD - h))
    return leg.prefactor / (roots.rD - h) *
        (_ellip_f(s, c, leg.m1) - alpha * _ellip_pole(s, c, leg.m1, n1))
end

function _d2_radial_residues(a, energy, lz)
    rp = _rplus(a)
    rm = _d2_rminus(a)
    pplus = energy * (rp^2 + a^2) - a * lz
    pminus = energy * (rm^2 + a^2) - a * lz
    return (
        rplus=rp,
        rminus=rm,
        c_phi_plus=a * pplus / (rp - rm),
        c_phi_minus=a * pminus / (rm - rp),
        c_t_plus=2 * rp * pplus / (rp - rm),
        c_t_minus=2 * rm * pminus / (rm - rp),
    )
end

function _d2_radial_phi_dot(a, energy, lz, r)
    return a * (energy * (r^2 + a^2) - a * lz) / (r^2 - 2 * r + a^2) -
           a * energy
end

function _d2_radial_time_dot(a, energy, lz, r)
    return (r^2 + a^2) * (energy * (r^2 + a^2) - a * lz) /
           (r^2 - 2 * r + a^2)
end

function _d1_roots(outcome)
    roots = collect(outcome.roots)
    length(roots) == 3 || return nothing
    sort!(roots)
    return (x1=roots[1], x2=roots[2], x3=roots[3])
end

function _d1_rminus(a)
    return 1 - sqrt(max(0.0, 1 - a^2))
end

function _d1_p(a, lz, r)
    return r^2 + a^2 - a * lz
end

function _d1_radial_residues(a, lz)
    rp = _rplus(a)
    rm = _d1_rminus(a)
    pplus = _d1_p(a, lz, rp)
    pminus = _d1_p(a, lz, rm)
    return (
        rplus=rp,
        rminus=rm,
        c_phi_plus=a * pplus / (rp - rm),
        c_phi_minus=a * pminus / (rm - rp),
        c_t_plus=2 * rp * pplus / (rp - rm),
        c_t_minus=2 * rm * pminus / (rm - rp),
    )
end

function _d1_positive_phi_radial_infinity(a, lz, leg)
    residues = _d1_radial_residues(a, lz)
    return residues.c_phi_plus * _three_real_pole(leg, residues.rplus, pi / 2) +
        residues.c_phi_minus * _three_real_pole(leg, residues.rminus, pi / 2)
end

function _d1_radial_phi_dot(a, lz, r)
    return a * _d1_p(a, lz, r) / (r^2 - 2 * r + a^2) - a
end

function _d1_radial_time_dot(a, lz, r)
    return (r^2 + a^2) * _d1_p(a, lz, r) / (r^2 - 2 * r + a^2)
end

function _unsupported_scatter(parameters, constants, outcome, reason)
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    return KerrGeoScatter(outcome.formula, outcome.energy_regime, outcome.outcome, parameters,
        constants, "Mino", (radial=outcome.roots, root_class=outcome.root_class),
        (kind=:closest_approach, lambda0=0.0, t0=0.0, phi0=0.0),
        nothing, nothing,
        (radial=r -> kerr_radial_potential(a, energy, lz, q, r),
         polar_z=z -> kerr_polar_z_potential(a, energy, lz, q, z)),
        (radial=nothing, polar_z=nothing),
        (supported=false, reason=reason))
end

"""
    kerr_geo_scatter(a, p, e, x; input=:apex, kwargs...)
    kerr_geo_scatter(a, E, Lz, Q; input=:constants, polar_phase=0.0)
    kerr_geo_scatter(a, constants; kwargs...)

Finite-window scatter orbit (E ≥ 1, from infinity to infinity through a turning
point). Inclined orbits need constants input; `polar_phase` is the polar argument at
closest approach. Constants without a scatter component return a record with
`Status.supported == false`.
"""
function kerr_geo_scatter(a::Real, b::Real, c::Real, d::Real; input::Symbol=:apex,
        constants=nothing, polar_phase=nothing)
    if input === :constants
        constants_tuple = (E=b, Lz=c, Q=d)
        parameters = (a=a, input=:constants)
    elseif input === :apex
        if constants === nothing
            c_ = kerr_geo_constants_of_motion(a, b, c, d)
            constants_tuple = (E=c_["E"], Lz=c_["Lz"], Q=c_["Q"])
        else
            constants_tuple = constants
        end
        parameters = (a=a, p=b, e=c, x=d, input=:apex)
    else
        error("Unknown scatter input type. Use :apex or :constants.")
    end
    outcome = _outcome_at_infinity(a, constants_tuple.E, constants_tuple.Lz, constants_tuple.Q)
    outcome.outcome === :scatter || return _unsupported_scatter(parameters,
        constants_tuple, outcome, "Outcome is $(outcome.outcome), not scatter.")
    input === :apex && !isapprox(maximum(outcome.roots), b / (1 + c); rtol=1e-6) &&
        return _unsupported_scatter(parameters, constants_tuple, outcome,
            "The periapsis p/(1+e) is not the turning point of the orbit from infinity.")
    return _scatter_finite_window(parameters, constants_tuple, outcome;
        polar_phase=polar_phase === nothing ? 0.0 : polar_phase)
end

function kerr_geo_scatter(a::Real, constants::NamedTuple; kwargs...)
    return kerr_geo_scatter(a, constants.E, constants.Lz, constants.Q; input=:constants, kwargs...)
end

function kerr_geo_scatter(a::Real, constants::Tuple{<:Real,<:Real,<:Real}; kwargs...)
    return kerr_geo_scatter(a, constants[1], constants[2], constants[3]; input=:constants, kwargs...)
end

function _impact_effective_magnitude(energy, lz, q)
    kerr_energy_regime(energy) === :hyperbolic || return Inf
    return sqrt(max(lz^2 + max(q, 0.0), 0.0)) / sqrt(_e2m1(energy))
end

function _clamp_unit(x)
    return clamp(x, -1.0, 1.0)
end

function _cartesian_er(q, phi)
    s = sqrt(max(1 - q^2, 0.0))
    return (x=s * cos(phi), y=s * sin(phi), z=q)
end

function _cartesian_etheta(q, phi)
    s = sqrt(max(1 - q^2, 0.0))
    return (x=q * cos(phi), y=q * sin(phi), z=-s)
end

function _cartesian_ephi(phi)
    return (x=-sin(phi), y=cos(phi), z=0.0)
end

function _tuple_scale(c, v)
    return (x=c * v.x, y=c * v.y, z=c * v.z)
end

function _tuple_dot(a, b)
    return a.x * b.x + a.y * b.y + a.z * b.z
end

function _tuple_norm(v)
    return sqrt(max(_tuple_dot(v, v), 0.0))
end

function _d2_infinity_amplitude(roots)
    return atan(sqrt((roots.rC - roots.rA) / (roots.rD - roots.rC)))
end

function _d2_positive_phi_radial_amplitude(a, energy, lz, leg, phi)
    residues = _d2_radial_residues(a, energy, lz)
    return residues.c_phi_plus * _d2_pole_primitive(leg, residues.rplus, phi) +
           residues.c_phi_minus * _d2_pole_primitive(leg, residues.rminus, phi)
end

function _d2_positive_total_phi_infinity(a, energy, lz, leg)
    phi_inf = _d2_infinity_amplitude(leg.roots)
    return _d2_positive_phi_radial_amplitude(a, energy, lz, leg, phi_inf) +
           lz * leg.lambda_infinity
end

function _d2_roots_from_scatter(kg::KerrGeoScatter)
    radial = kg.Roots.radial
    length(radial) == 4 || error("D2 diagnostics require four radial roots.")
    return (rA=radial[1], rB=radial[2], rC=radial[3], rD=radial[4])
end

function _d2_equatorial_asymptotic_diagnostics(kg::KerrGeoScatter)
    a = kg.OrbitalParameters.a
    energy = kg.ConstantsOfMotion.E
    lz = kg.ConstantsOfMotion.Lz
    roots = _d2_roots_from_scatter(kg)
    p_inf = sqrt(_e2m1(energy))
    impact_magnitude = abs(lz) / p_inf
    phi_turn_to_infinity = _d2_positive_total_phi_infinity(a, energy, lz, _four_real_leg(energy, roots))
    delta_phi = 2 * phi_turn_to_infinity
    lz_sign = sign(iszero(lz) ? 1.0 : lz)
    signed_deflection = delta_phi - lz_sign * pi
    unsigned_deflection = abs(delta_phi) - pi
    principal_cosine = -cos(delta_phi)
    principal_angle = acos(_clamp_unit(principal_cosine))

    phi_in = -phi_turn_to_infinity
    phi_out = phi_turn_to_infinity
    er_in = _cartesian_er(0.0, phi_in)
    er_out = _cartesian_er(0.0, phi_out)
    n_in = _tuple_scale(-1.0, er_in)
    n_out = er_out

    alpha = -lz / p_inf
    beta = 0.0
    ephi_in = _cartesian_ephi(phi_in)
    etheta_in = _cartesian_etheta(0.0, phi_in)
    impact_vector_raw = (
        x=alpha * ephi_in.x + beta * etheta_in.x,
        y=alpha * ephi_in.y + beta * etheta_in.y,
        z=alpha * ephi_in.z + beta * etheta_in.z,
    )
    impact_norm = _tuple_norm(impact_vector_raw)

    return (
        formula=kg.Formula,
        supported_trajectory=kg.Status.supported,
        impact_vector=(x=impact_vector_raw.x,
                       y=impact_vector_raw.y,
                       z=impact_vector_raw.z,
                       status=:equatorial_hyperbolic_screen,
                       frame=:incoming_asymptotic_screen,
                       dot_incoming_direction=_tuple_dot(impact_vector_raw, n_in)),
        impact_magnitude=(value=impact_norm,
                          status=:equatorial_hyperbolic_scalar,
                          expected=impact_magnitude,
                          screen_coordinates=(alpha=alpha, beta=beta)),
        incoming_direction=(value=(x=n_in.x, y=n_in.y, z=n_in.z),
                            status=:equatorial_hyperbolic_direction,
                            norm=_tuple_norm(n_in),
                            position_phi=phi_in),
        outgoing_direction=(value=(x=n_out.x, y=n_out.y, z=n_out.z),
                            status=:equatorial_hyperbolic_direction,
                            norm=_tuple_norm(n_out),
                            position_phi=phi_out),
        azimuthal_deflection=(value=signed_deflection,
                              unsigned=unsigned_deflection,
                              total_azimuth=delta_phi,
                              turn_to_infinity=phi_turn_to_infinity,
                              status=:equatorial_hyperbolic_unwrapped,
                              sign_convention=:delta_phi_minus_sign_lz_pi),
        deflection_angle_3d=(value=principal_angle,
                             cosine=principal_cosine,
                             status=:equatorial_hyperbolic_principal),
        convention=(
            straight_line_zero=:momentum_direction,
            position_ray_angle_is_distinct=true,
            impact_vector_convention=:incoming_screen_alpha_ephi_plus_beta_etheta,
            phi_reference=:closest_approach_zero,
        ),
    )
end

function _d2_generic_hyperbolic_asymptotic_diagnostics(kg::KerrGeoScatter)
    a = kg.OrbitalParameters.a
    energy = kg.ConstantsOfMotion.E
    lz = kg.ConstantsOfMotion.Lz
    qcarter = kg.ConstantsOfMotion.Q
    roots = _d2_roots_from_scatter(kg)
    lambda_infinity = kg.Status.domain.mino[2]
    p_inf = sqrt(_e2m1(energy))

    polar_phase = kg.ReferenceZero.polar_phase === nothing ? 0.0 : kg.ReferenceZero.polar_phase
    # the trajectory's own polar motion (same phase convention), φ zero at closest approach
    polar = _window_polar_motion(a, energy, lz, qcarter, polar_phase)
    polar_in = (z=polar.z(-lambda_infinity), uz=polar.uz(-lambda_infinity),
        phi=polar.phi(-lambda_infinity))
    polar_out = (z=polar.z(lambda_infinity), uz=polar.uz(lambda_infinity),
        phi=polar.phi(lambda_infinity))
    radial_phi_infinity = _d2_positive_phi_radial_amplitude(
        a, energy, lz, _four_real_leg(energy, roots), _d2_infinity_amplitude(roots))
    phi_in = -radial_phi_infinity + polar_in.phi
    phi_out = radial_phi_infinity + polar_out.phi
    delta_phi = phi_out - phi_in

    z_in = polar_in.z
    z_out = polar_out.z
    sin_in = sqrt(max(1 - z_in^2, 0.0))
    sin_in > 0 || error("Generic D2 asymptotic impact screen requires sin(theta_in)>0.")
    theta_dot_in = -polar_in.uz / sin_in
    alpha = -lz / (p_inf * sin_in)
    beta = theta_dot_in / p_inf

    er_in = _cartesian_er(z_in, phi_in)
    er_out = _cartesian_er(z_out, phi_out)
    n_in = _tuple_scale(-1.0, er_in)
    n_out = er_out
    principal_cosine = _tuple_dot(n_in, n_out)
    principal_angle = acos(_clamp_unit(principal_cosine))

    ephi_in = _cartesian_ephi(phi_in)
    etheta_in = _cartesian_etheta(z_in, phi_in)
    impact_vector_raw = (
        x=alpha * ephi_in.x + beta * etheta_in.x,
        y=alpha * ephi_in.y + beta * etheta_in.y,
        z=alpha * ephi_in.z + beta * etheta_in.z,
    )
    impact_norm = _tuple_norm(impact_vector_raw)
    theta_potential_in = kerr_polar_z_potential(a, energy, lz, qcarter, z_in) / (sin_in^2)
    impact_identity = sqrt(max(
        (qcarter + lz^2 + a^2 * _e2m1(energy) * z_in^2) / _e2m1(energy),
        0.0,
    ))

    return (
        formula=kg.Formula,
        supported_trajectory=kg.Status.supported,
        impact_vector=(x=impact_vector_raw.x,
                       y=impact_vector_raw.y,
                       z=impact_vector_raw.z,
                       status=:generic_hyperbolic_screen,
                       frame=:incoming_asymptotic_screen,
                       dot_incoming_direction=_tuple_dot(impact_vector_raw, n_in)),
        impact_magnitude=(value=impact_norm,
                          status=:generic_hyperbolic_screen_norm,
                          expected=impact_identity,
                          screen_coordinates=(alpha=alpha, beta=beta),
                          theta_potential=theta_potential_in),
        incoming_direction=(value=(x=n_in.x, y=n_in.y, z=n_in.z),
                            status=:generic_hyperbolic_direction,
                            norm=_tuple_norm(n_in),
                            position_phi=phi_in,
                            q=z_in),
        outgoing_direction=(value=(x=n_out.x, y=n_out.y, z=n_out.z),
                            status=:generic_hyperbolic_direction,
                            norm=_tuple_norm(n_out),
                            position_phi=phi_out,
                            q=z_out),
        azimuthal_deflection=(value=delta_phi,
                              total_azimuth=delta_phi,
                              status=:generic_hyperbolic_unwrapped_phase,
                              sign_convention=:unwrapped_bl_phase_not_scalar_deflection),
        deflection_angle_3d=(value=principal_angle,
                             cosine=principal_cosine,
                             status=:generic_hyperbolic_principal),
        convention=(
            straight_line_zero=:momentum_direction,
            position_ray_angle_is_distinct=true,
            impact_vector_convention=:incoming_screen_alpha_ephi_plus_beta_etheta,
            phi_reference=:closest_approach_zero,
            polar_phase=kg.ReferenceZero.polar_phase,
            polar_phase_convention=polar.metadata.phase_convention,
        ),
    )
end

function _d1_roots_from_scatter(kg::KerrGeoScatter)
    radial = kg.Roots.radial
    length(radial) == 3 || error("D1 diagnostics require three radial roots.")
    return (x1=radial[1], x2=radial[2], x3=radial[3])
end

function _d1_parabolic_geometric_diagnostics(kg::KerrGeoScatter)
    a = kg.OrbitalParameters.a
    lz = kg.ConstantsOfMotion.Lz
    qcarter = kg.ConstantsOfMotion.Q
    leg = _three_real_leg(_d1_roots_from_scatter(kg))
    lambda_infinity = leg.lambda_infinity
    polar_phase = kg.ReferenceZero.polar_phase === nothing ? 0.0 : kg.ReferenceZero.polar_phase
    generic_inclined = !iszero(qcarter)

    radial_phi_infinity = _d1_positive_phi_radial_infinity(a, lz, leg)
    if generic_inclined
        # the trajectory's own polar motion (same phase convention)
        polar = _window_polar_motion(a, kg.ConstantsOfMotion.E, lz, qcarter, polar_phase)
        polar_in = (z=polar.z(-lambda_infinity), phi=polar.phi(-lambda_infinity))
        polar_out = (z=polar.z(lambda_infinity), phi=polar.phi(lambda_infinity))
        z_in = polar_in.z
        z_out = polar_out.z
        phi_in = -radial_phi_infinity + polar_in.phi
        phi_out = radial_phi_infinity + polar_out.phi
    else
        z_in = 0.0
        z_out = 0.0
        phi_turn_to_infinity = radial_phi_infinity + lz * lambda_infinity
        phi_in = -phi_turn_to_infinity
        phi_out = phi_turn_to_infinity
    end

    delta_phi = phi_out - phi_in
    er_in = _cartesian_er(z_in, phi_in)
    er_out = _cartesian_er(z_out, phi_out)
    n_in = _tuple_scale(-1.0, er_in)
    n_out = er_out
    principal_cosine = _tuple_dot(n_in, n_out)
    principal_angle = acos(_clamp_unit(principal_cosine))
    angular_scale = sqrt(max(lz^2 + max(qcarter, 0.0), 0.0))
    lz_sign = sign(iszero(lz) ? 1.0 : lz)
    signed_equatorial = generic_inclined ? NaN : delta_phi - lz_sign * pi
    unsigned_equatorial = generic_inclined ? NaN : abs(delta_phi) - pi

    return (
        formula=kg.Formula,
        supported_trajectory=kg.Status.supported,
        impact_vector=(x=NaN,
                       y=NaN,
                       z=NaN,
                       status=:parabolic_hyperbolic_screen_singular,
                       reason=:p_infinity_zero),
        impact_magnitude=(value=angular_scale,
                          status=:parabolic_angular_scale,
                          hyperbolic_impact=:singular,
                          strict_infinity_impact=:not_defined),
        incoming_direction=(value=(x=n_in.x, y=n_in.y, z=n_in.z),
                            status=:parabolic_geometric_direction,
                            norm=_tuple_norm(n_in),
                            position_phi=phi_in,
                            q=z_in),
        outgoing_direction=(value=(x=n_out.x, y=n_out.y, z=n_out.z),
                            status=:parabolic_geometric_direction,
                            norm=_tuple_norm(n_out),
                            position_phi=phi_out,
                            q=z_out),
        azimuthal_deflection=(value=generic_inclined ? delta_phi : signed_equatorial,
                              unsigned=unsigned_equatorial,
                              total_azimuth=delta_phi,
                              status=generic_inclined ?
                                     :generic_parabolic_unwrapped_phase :
                                     :equatorial_parabolic_unwrapped_geometric,
                              sign_convention=generic_inclined ?
                                              :unwrapped_bl_phase_not_scalar_deflection :
                                              :delta_phi_minus_sign_lz_pi),
        deflection_angle_3d=(value=principal_angle,
                             cosine=principal_cosine,
                             status=:parabolic_geometric_principal),
        convention=(
            straight_line_zero=:geometric_asymptote_direction,
            position_ray_angle_is_distinct=true,
            impact_vector_convention=:not_available_for_parabolic_p_infinity_zero,
            phi_reference=:closest_approach_zero,
            polar_phase=kg.ReferenceZero.polar_phase,
            parabolic_scale=:sqrt_lz2_plus_q,
        ),
    )
end

"""
    kerr_geo_scatter_asymptotic_diagnostics(orbit)

Asymptotic data of a scatter orbit, with φ = 0 at closest approach: the incoming and outgoing
directions, the azimuthal deflection and the 3D angle between the asymptotic directions. D2
(E > 1) also gives the impact vector on the incoming screen and its magnitude. The azimuthal
deflection is Δφ − sign(Lz)π for equatorial orbits and the unwrapped Boyer–Lindquist phase Δφ
for inclined ones. For D1 (E = 1) the directions are the geometric asymptote directions, the
impact vector is undefined because √(E² − 1) = 0, and `impact_magnitude` is the angular scale
√(Lz² + Q).
"""
function kerr_geo_scatter_asymptotic_diagnostics(kg::KerrGeoScatter)
    kg.Outcome === :scatter || error("Asymptotic diagnostics require a scatter orbit.")
    energy = kg.ConstantsOfMotion.E
    lz = kg.ConstantsOfMotion.Lz
    q = kg.ConstantsOfMotion.Q
    equatorial = iszero(q)
    regime = kerr_energy_regime(energy)
    hyperbolic = regime === :hyperbolic
    if kg.Formula === :hyperbolic_scatter && equatorial && hyperbolic && kg.Status.supported
        return _d2_equatorial_asymptotic_diagnostics(kg)
    elseif kg.Formula === :hyperbolic_scatter && !equatorial && hyperbolic && kg.Status.supported
        return _d2_generic_hyperbolic_asymptotic_diagnostics(kg)
    elseif kg.Formula === :parabolic_scatter && regime === :parabolic && kg.Status.supported
        return _d1_parabolic_geometric_diagnostics(kg)
    end
    magnitude = _impact_effective_magnitude(energy, lz, q)
    vector = if equatorial && isfinite(magnitude)
        (x=0.0, y=sign(lz == 0 ? 1.0 : lz) * magnitude, z=0.0,
         status=:equatorial_hyperbolic_screen)
    elseif equatorial
        (x=0.0, y=sign(lz == 0 ? 1.0 : lz) * Inf, z=0.0,
         status=:parabolic_velocity_at_infinity_singular)
    else
        (x=NaN, y=NaN, z=NaN,
         status=:unavailable)
    end
    scalar_status =
        hyperbolic && equatorial ? :equatorial_hyperbolic_scalar :
        hyperbolic ? :effective_impact_norm :
        :parabolic_velocity_at_infinity_singular
    return (
        formula=kg.Formula,
        supported_trajectory=kg.Status.supported,
        impact_vector=vector,
        impact_magnitude=(value=magnitude, status=scalar_status),
        incoming_direction=(value=nothing, status=:unavailable),
        outgoing_direction=(value=nothing, status=:unavailable),
        azimuthal_deflection=(value=NaN, status=:unavailable),
        deflection_angle_3d=(value=NaN, status=:unavailable),
        convention=(
            straight_line_zero=:momentum_direction,
            position_ray_angle_is_distinct=true,
        ),
    )
end

"""
    kerr_geo_scatter_asymptotic_state(orbit, side)

Asymptotic endpoint of a scatter orbit on `side` (`:incoming` = past infinity, `:outgoing` =
future infinity): its Mino-time endpoint and the data of `kerr_geo_scatter_asymptotic_diagnostics`.
"""
function kerr_geo_scatter_asymptotic_state(kg::KerrGeoScatter, side::Symbol)
    side in (:incoming, :outgoing) ||
        error("Asymptotic state side must be :incoming or :outgoing.")
    diagnostics = kerr_geo_scatter_asymptotic_diagnostics(kg)
    return (
        formula=kg.Formula,
        side=side,
        endpoint_role=side === :incoming ? :past_infinity : :future_infinity,
        ordinary_trajectory_sample=false,
        mino_endpoint=kg.Status.domain.mino[side === :incoming ? 1 : 2],
        diagnostics=diagnostics,
        convention=diagnostics.convention,
    )
end

function Base.show(io::IO, kg::KerrGeoScatter)
    print(io, "KerrGeoScatter(", kg.Formula, ", constants=")
    show(io, kg.ConstantsOfMotion)
    print(io, ", supported=", kg.Status.supported, ")")
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoScatter)
    println(io, "KerrGeoScatter(")
    print(io, "    Formula = "); show(io, kg.Formula); println(io, ",")
    print(io, "    EnergyRegime = "); show(io, kg.EnergyRegime); println(io, ",")
    print(io, "    Outcome = "); show(io, kg.Outcome); println(io, ",")
    print(io, "    ConstantsOfMotion = "); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    ReferenceZero = "); show(io, kg.ReferenceZero); println(io, ",")
    print(io, "    Status = "); show(io, kg.Status); println(io, ",")
    print(io, ")")
end

# Radial closed forms (r(λ), λ(r), rates) of the two finite-window scatter formulas: D2
# (E > 1, four real roots) and D1 (E = 1, three real roots). λ = 0 at the turning point.
function _scatter_radial_model(formula, a, energy, lz, outcome)
    if formula === :hyperbolic_scatter
        roots = _d2_roots(outcome)
        roots === nothing && return nothing
        leg = _four_real_leg(energy, roots)
        λ∞ = leg.lambda_infinity
        return (radius=s -> (t = _leg_time(s, λ∞);
                _four_real_radius(leg, t / leg.prefactor, (λ∞ - t) / leg.prefactor)),
            tdot=r -> _d2_radial_time_dot(a, energy, lz, r),
            phidot=r -> _d2_radial_phi_dot(a, energy, lz, r),
            mino=r -> _four_real_mino_from_turn(leg, r),
            lambda_infinity=λ∞, r_turn=roots.rD,
            radial=(roots.rA, roots.rB, roots.rC, roots.rD), modulus=leg.m)
    else
        roots = _d1_roots(outcome)
        roots === nothing && return nothing
        leg = _three_real_leg(roots)
        λ∞ = leg.lambda_infinity
        return (radius=s -> _three_real_radius_from_infinity(leg, λ∞ - _leg_time(s, λ∞)),
            tdot=r -> _d1_radial_time_dot(a, lz, r),
            phidot=r -> _d1_radial_phi_dot(a, lz, r),
            mino=r -> λ∞ - _three_real_mino_from_infinity(leg, r),
            lambda_infinity=λ∞, r_turn=roots.x3,
            radial=(roots.x1, roots.x2, roots.x3), modulus=leg.m)
    end
end

"""
Finite-window scatter orbit for E > 1 (four real roots) or E = 1 (three real roots),
symmetric about the turning point at λ = 0 (λ < 0 incoming, λ > 0 outgoing). Equatorial
orbits and generic inclined ones with a fixed `polar_phase` at closest approach are
supported.
"""
function _scatter_finite_window(parameters, constants, outcome; polar_phase=0.0)
    formula = outcome.formula
    a, energy, lz, q = parameters.a, constants.E, constants.Lz, constants.Q
    # The companion roots need the same metric-based correction as component roots.
    coefficients = kerr_radial_coefficients(a, energy, lz, q)
    corrected_roots = sort([_polish_root(coefficients, r) for r in outcome.roots])
    outcome = merge(outcome, (roots=Tuple(corrected_roots),))
    unsupported(msg) = _unsupported_scatter(parameters, constants, outcome, msg)
    radial = _scatter_radial_model(formula, a, energy, lz, outcome)
    radial === nothing && return unsupported(
        "The radial roots do not have the structure required by the $formula formula.")
    polar = _window_polar_motion(a, energy, lz, q, polar_phase)
    inclined = polar.inclined
    inclined && parameters.input !== :constants && return unsupported(
        "Inclined scatter orbits need constants input (the APEX polar-phase convention is not defined).")

    R(r) = kerr_radial_potential(a, energy, lz, q, r)
    Θ(z) = kerr_polar_z_potential(a, energy, lz, q, z)
    r_of_lambda(λ) = radial.radius(abs(λ))
    branch(λ) = λ < 0 ? -1.0 : 1.0
    rdot(λ) = branch(λ) * sqrt(max(R(r_of_lambda(λ)), 0.0))
    # t, φ: radial spectral engine (infinity → turning point → infinity) + polar primitive,
    # zero at the turning point
    coords = _engine_coordinates(a, energy, lz, q, r_of_lambda, polar.primitive;
        potential=_coefficient_potential(a, energy, lz, q),
        domain=(-radial.lambda_infinity, radial.lambda_infinity), ends=(:infinity, :infinity),
        turn=0.0, σ=1.0, λ_bl=0.0)
    t(λ) = _coords_t(coords, λ)
    phi(λ) = _coords_phi(coords, λ)
    tdot(λ) = radial.tdot(r_of_lambda(λ)) + polar.tdot(λ)
    phidot(λ) = radial.phidot(r_of_lambda(λ)) + polar.phidot(λ)
    utheta(λ) = inclined ? -polar.uz(λ) / sqrt(max(1 - polar.z(λ)^2, 0.0)) : 0.0
    function sample_by_radius(r; branch=:outgoing)
        λ = radial.mino(r)
        branch === :incoming && return (lambda=-λ, r=r)
        branch === :outgoing && return (lambda=λ, r=r)
        error("sample_by_radius branch must be :incoming or :outgoing.")
    end
    phase = _polar_phase_metadata(polar, polar_phase, :closest_approach)

    return KerrGeoScatter(formula, outcome.energy_regime, :scatter, parameters, constants, "Mino",
        (radial=radial.radial, root_class=outcome.root_class, modulus=radial.modulus,
         polar=inclined ? polar.rootdata : nothing),
        (kind=:closest_approach, lambda0=0.0, t0=0.0, phi0=0.0, r_turn=radial.r_turn,
         polar_phase=inclined ? polar_phase : nothing),
        (t=t, r=r_of_lambda, theta=polar.theta, phi=phi, q=polar.z,
         sample_by_radius=sample_by_radius),
        (ut=tdot, ur=rdot, utheta=utheta, uphi=phidot),
        (radial=R, polar_z=Θ),
        (radial=λ -> rdot(λ)^2 - R(r_of_lambda(λ)),
         polar_z=λ -> polar.uz(λ)^2 - Θ(polar.z(λ))),
        (supported=true, reason="ok", domain=_scatter_domain(radial.lambda_infinity),
         polar_phase=phase))
end
