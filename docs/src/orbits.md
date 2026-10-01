# Working with an orbit

## The family

[`kerr_geodesic`](@ref) takes the spin and the constants of motion, as a tuple or a
NamedTuple, or the spin and the APEX parameters ``(p, e, x)`` of an orbit with a periapsis:

```julia
kerr_geodesic(a, (E, Lz, Q))
kerr_geodesic(a, (E = E, Lz = Lz, Q = Q))
kerr_geodesic(a, p, e, x)
```

The result is a [`KerrGeodesicFamily`](@ref) that holds every orbit the constants allow:

| Field | Content |
| :--- | :--- |
| `Stable`, `Plunge`, `Capture`, `Scatter`, `Trapped` | the member of that class, or `nothing` |
| `Critical` | a tuple of Critical members, ordered by role (on the root, outer, inner) |
| `ConstantsOfMotion` | `(E, Lz, Q)` |
| `InputType`, `Parameters` | `:constants` with `(a,)`, or `:apex` with `(a, p, e, x)` |
| `BroadClass` | one class that summarises the family, `:critical` whenever it has Critical members |
| `RootClass` | a label for the root structure of ``R`` |
| `Status` | `supported`, `reason`, `case_ids` and the classification details |

Default display shows a short summary of the family or member, not its full diagnostic
records or function objects. The fields remain available, for example `kg.Status.reason`
and `m.Status.spectral.achieved`. Display does not evaluate the trajectory or build
spectral tables.

```@example orbits
using KerrGeodesics

kg = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
kg.Status.case_ids
```

Constants that allow no motion outside the horizon give a family with no members,
`Status.supported == false` and a `Status.reason`.

## The member

Every member is a [`KerrGeoComponent`](@ref)`{C}`, with its class `C` as a type parameter,
and has the same fields in every class:

| Field | Content |
| :--- | :--- |
| `CaseId` | the case, such as `:A1` or `:K4`, or the tier member, such as `:A_H1` |
| `Tier` | `:primary`, `:horizon` or `:extremal` |
| `Role` | `:on_root`, `:outer` or `:inner` for Critical members, `:none` otherwise |
| `ConstantsOfMotion` | `(a, E, Lz, Q)` |
| `Roots` | `radial`, root data whose representation depends on the member: real-root tuples, complex-pair parameters or records with multiplicity; `polar`, the polar solution (its sector, the turning values of ``z^2``, the phase convention) |
| `Domain` | `mino`, the range of ``λ``, with the role of each end; `horizon_lambda` for members that end on the future horizon |
| `ReferenceZero` | the events where ``λ``, ``t``, ``φ``, ``τ``, ``v``, ``ψ`` vanish (see [Where the coordinates are zero](@ref)) |
| `Trajectory` | the coordinates, as functions of ``λ`` |
| `Velocity` | the rates ``dx^μ/dλ`` and ``dτ/dλ``, as functions of ``λ`` |
| `Potentials` | the potentials as functions of their own variable: `radial(r)` ``= R(r)``, `polar_z(z)` ``= Θ(z)`` |
| `Residuals` | functions of ``λ``: `radial` ``= (dr/dλ)^2 - R``, `polar_z` ``= (dz/dλ)^2 - Θ``, `normalization` ``= g_{μν}u^μu^ν + 1`` |
| `Status` | facts about the member: `supported`, `spectral` (see [Numerics and accuracy](@ref)), and for Stable members `apex`, `frequencies` and `precision` |
| `Component` | the classified region as a [`KerrGeoRadialComponent`](@ref), or `nothing` for members built without one |

`Roots.radial` is not a uniform array of radii. For example, a Trapped N4 member gives
real roots `x1,x2`, a complex-pair description `rho,eta`, and auxiliary distances `A,B`.
Exact-extremal primary members may give records with `radius`, `multiplicity` and
residual metadata. Inspect the member's root data before treating its entries as radii.

The class constructors, such as [`kerr_geo_plunge_component`](@ref) or
[`kerr_geo_trapped`](@ref), build one member directly from the constants.

## Trajectory

`Trajectory` holds ``t``, ``r``, ``θ`` (`theta`), ``z = \cos θ``, ``φ`` (`phi`) and ``τ``
(`tau`) as functions of ``λ``:

```@example orbits
stable = kg.Stable
λ = 2.5
(t = stable.Trajectory.t(λ), r = stable.Trajectory.r(λ), theta = stable.Trajectory.theta(λ),
 z = stable.Trajectory.z(λ), phi = stable.Trajectory.phi(λ), tau = stable.Trajectory.tau(λ))
```

They accept ``λ`` in `Domain.mino`, except at endpoints where that coordinate diverges.
A stable orbit runs forever. In the following subextremal example the plunge starts at
its turning point at ``λ = 0``; axis infall and exact-extremal crossing use other origins
(see [Where the coordinates are zero](@ref)), and exact-critical plunges approach the
horizon only asymptotically:

```@example orbits
stable.Domain.mino, kg.Plunge.Domain.mino
```

Outside this range the functions throw a `DomainError`, and so do ``t`` and ``φ`` at an end
on a horizon, where they diverge. Besides these six, a member carries the functions its
motion calls for, and `keys(member.Trajectory)` lists them:

| Functions | What they give | Typical members |
| :--- | :--- | :--- |
| `v`, `psi` | the ingoing coordinates ``v``, ``ψ``, finite on the future horizon | members that cross the future horizon, K7 and K10, and members at ``a = ±1`` |
| `rstar` | the tortoise coordinate ``r_*`` | the same, except Trapped members |
| `u`, `chi` | the outgoing coordinates ``u``, ``χ``, finite on the past horizon | Trapped members and members at ``a = ±1`` |
| `lambda_of_radius` | ``λ`` at a given radius | Critical members off the root, and Plunge, Capture, Scatter and Trapped members for ``a^2 < 1`` |
| `radial_mino_increment`, `radial_time_increment`, `radial_phi_increment`, `radial_proper_increment` | radial integrals between two radii | the same, except Trapped members |
| `radial_t`, `radial_phi`, `radial_tau` | the radial parts of ``t``, ``φ``, ``τ`` | the same, except Trapped members |
| `radial_v_increment`, `radial_psi_increment` | radial integrals for ``v`` and ``ψ`` | K2, K5, K7, K8, K10, K11 and Capture members |
| `full`, `outgoing`, `incoming` | all coordinates at one ``λ``, as a NamedTuple | Trapped members |

## Four-velocity

`Velocity` holds the Mino-time rates `ut` ``= dt/dλ``, `ur` ``= dr/dλ``, `uz` ``= dz/dλ``,
`utheta` ``= dθ/dλ``, `uphi` ``= dφ/dλ``, and `dtau_dlambda` ``= Σ``. The four-velocity
``u^μ = dx^μ/dτ`` is their ratio:

```@example orbits
rates = stable.Velocity
four_velocity = (rates.ut(λ), rates.ur(λ), rates.utheta(λ), rates.uphi(λ)) ./ rates.dtau_dlambda(λ)
```

`Residuals.normalization` is ``g_{μν}u^μu^ν + 1``, zero for a timelike geodesic up to
rounding:

```@example orbits
stable.Residuals.normalization(λ)
```

For a bound orbit given by ``(a, p, e, x)``, [`kerr_geo_four_velocity`](@ref) returns the
four components as functions of ``λ`` directly, contravariant or, with `Covariant = true`,
covariant.

## Orbital parameters

```@example orbits
stable.ConstantsOfMotion
```

A Stable member also knows its APEX parameters, computed from its turning points:
``r = p/(1 - e)`` and ``r = p/(1 + e)`` for the radial motion and
``z_{max} = \sqrt{1 - x^2}`` for the polar motion.

```@example orbits
stable.Status.apex
```

```@example orbits
stable.Roots.radial
```

```@example orbits
sqrt(stable.Roots.polar.zminus)
```

## Frequencies

A Stable member gives its Mino-time frequencies ``Υ_r``, ``Υ_θ``, ``Υ_φ`` and the mean
rate ``Υ_t = ⟨dt/dλ⟩``:

```@example orbits
ϒ = stable.Status.frequencies
```

The Boyer–Lindquist frequencies, the frequencies in coordinate time ``t`` seen from infinity,
are ``Ω_i = Υ_i / Υ_t``:

```@example orbits
(Ωr = ϒ.ϒr / ϒ.ϒt, Ωθ = ϒ.ϒθ / ϒ.ϒt, Ωφ = ϒ.ϒϕ / ϒ.ϒt)
```

[`kerr_geo_frequencies`](@ref) computes the same from ``(a, p, e, x)``, in Mino time,
Boyer–Lindquist time or proper time:

```@example orbits
kerr_geo_frequencies(0.9, 10.0, 0.5, 0.8; Time = "BoyerLindquist")
```

## At the horizon

Boyer–Lindquist ``t`` and ``φ`` diverge where an orbit crosses the horizon. Members that
cross it also carry the tortoise coordinate ``r_*`` and the ingoing coordinates
``v = t + r_*`` and ``ψ = φ + φ_H``. The latter two, not ``r_*``, stay finite there and
use the zeros specified by `ReferenceZero` (see [Coordinates regular at the horizon](@ref)).

```@example orbits
plunge = kerr_geodesic(0.9, (0.94, 0.1, 12.0)).Plunge
λH = plunge.Domain.horizon_lambda
(r = plunge.Trajectory.r(λH), rplus = kerr_horizons(0.9).rplus,
 v = plunge.Trajectory.v(λH), psi = plunge.Trajectory.psi(λH))
```

Approaching ``λ_H``, ``t`` grows without bound while ``v`` settles to zero:

```@example orbits
[(δ, plunge.Trajectory.t(λH - δ), plunge.Trajectory.v(λH - δ)) for δ in (1e-2, 1e-5, 1e-8)]
```

## Radius as the variable

Members whose radius changes monotonically on each leg invert ``r(λ)``:

```@example orbits
λ₂ = plunge.Trajectory.lambda_of_radius(2.5)
(λ₂, plunge.Trajectory.r(λ₂))
```

The radial integrals between two radii,

```math
\int_{r_1}^{r_2} \frac{f(r)\,dr}{\sqrt{R(r)}}, \qquad
f = 1,\ \frac{(r^2+a^2)P(r)}{Δ},\ \frac{aP(r)}{Δ} - aE,\ r^2,
```

are `radial_mino_increment`, `radial_time_increment`, `radial_phi_increment` and
`radial_proper_increment`: the radial parts of the changes in ``λ``, ``t``, ``φ`` and ``τ``.
Along a leg on which ``r`` decreases, the change along the orbit is minus the integral.

```@example orbits
(plunge.Trajectory.radial_mino_increment(3.0, 2.0),
 plunge.Trajectory.lambda_of_radius(2.0) - plunge.Trajectory.lambda_of_radius(3.0))
```

These functions take the radii themselves, so they work at radii far beyond those ``λ`` can
resolve (see [Numerics and accuracy](@ref)).

## Many points at once

[`kerr_geo_sample`](@ref) evaluates the trajectory and the rates at many values of ``λ`` in
one call, faster than calling each function in a loop:

```@example orbits
sample = kerr_geo_sample(stable, range(0, 20; length = 5))
sample.r
```

The result is a NamedTuple of vectors: `lambda`, `t`, `r`, `theta`, `phi`, `tau`, `ut`,
`ur`, `utheta`, `uphi`.

## Keywords

`kerr_geodesic` and the class constructors take these keywords.

| Keyword | Effect |
| :--- | :--- |
| `polar_sector` | the polar sector, where the constants allow more than one |
| `polar_hemisphere` | `:north` (default) or `:south`, for polar motion confined to one hemisphere |
| `polar_phase` | the polar phase at the reference event, default 0 (see [The polar phase](@ref)) |
| `axis`, `phi0` | motion along the spin axis, `:north` or `:south`, at azimuth `phi0` |
| `reference_radius` | the radius where ``t`` and ``φ`` vanish on Critical and Capture members |
| `initPhases` | `(qt0, qr0, qθ0, qφ0)`, the phases of the Stable member at ``λ = 0`` |
| `trapped_component` | `:full` (default), `:outgoing` or `:incoming` part of the Trapped member |
| `case_id`, `initial_radius`, `radial_sign`, `endpoint_intent` | require that exactly one member matches, and record it in `Status.selected_case` |

`radial_sign` is `:inward`, `:outward` or `:zero` (constant radius); `endpoint_intent` is a
class Symbol, `:horizon`, `:infinity` or `:finite`. These four keywords keep every member in
the family; a selection that matches no member, or more than one, is an error.
