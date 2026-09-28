# APEX and finite-window interfaces

Besides `kerr_geodesic`, the package has two groups of functions with their own inputs and
record types. The APEX interface works with the parameters ``(a, p, e, x)`` of bound orbits
and returns the classic closed-form solution in Mino time. The finite-window constructors
build one plunge, capture or scattered orbit from the constants. Both describe the same
geodesics as the members of `kerr_geodesic`.

```@setup apex
using KerrGeodesics
```

## APEX parameters

[`kerr_geo_constants_of_motion`](@ref) converts ``(a, p, e, x)`` to the constants,
[`kerr_geo_radial_roots`](@ref) and [`kerr_geo_polar_roots`](@ref) give the roots of the
radial and polar potentials, and [`kerr_geo_frequencies`](@ref) the fundamental
frequencies in Mino time, Boyer–Lindquist time or proper time:

```@example apex
kerr_geo_constants_of_motion(0.9, 10.0, 0.5, 0.8)
```

```@example apex
kerr_geo_frequencies(0.9, 10.0, 0.5, 0.8; Time = "Mino")
```

The landmarks of bound motion are functions of the inclination: the innermost stable
spherical orbit [`kerr_geo_isso`](@ref) (the ISCO [`kerr_geo_isco`](@ref) for
equatorial orbits), the innermost bound spherical orbit [`kerr_geo_ibso`](@ref), and the
separatrix [`kerr_geo_separatrix`](@ref), below which an orbit of given ``e`` and ``x``
plunges:

```@example apex
(isco = kerr_geo_isco(0.9, 1.0), isso = kerr_geo_isso(0.9, 0.8),
 ibso = kerr_geo_ibso(0.9, 0.8), separatrix = kerr_geo_separatrix(0.9, 0.5, 0.8))
```

[`kerr_geo_orbit_type_metadata`](@ref) classifies an APEX point without building it:

```@example apex
meta = kerr_geo_orbit_type_metadata(0.9, 10.0, 0.5, 0.8)
(meta.family, meta.labels, meta.stability, meta.energy_regime)
```

## The closed-form bound orbit

[`kerr_geo_orbit`](@ref) returns the stable orbit, or the constant-radius Critical orbit with
``E < 1``, as a `Dict` in the form of Fujita and Hikida: the trajectory, the four-velocity,
the frequencies, the oscillating parts of ``t`` and ``φ`` (`"CrossFunction"`) and the
constants.

```@example apex
orbit = kerr_geo_orbit(0.9, 10.0, 0.5, 0.8)
t, r, θ, φ = orbit["Trajectory"]
(r(0.0), orbit["Frequencies"]["ϒr"], orbit["Stability"])
```

On a constant-radius Critical orbit, the ISCO, the ISSO or an unstable circular or spherical
orbit with ``E < 1``, the radial frequency is zero, and an unstable orbit also carries
`"RadialLyapunovExponent"`, the Mino-time growth rate of a radial perturbation. Every other
APEX point is an `ArgumentError` that names the class; `kerr_geodesic(a, p, e, x)` builds
it. [`kerr_geo_stable`](@ref) returns the same stable orbit as a `KerrGeoStable` record, and
[`kerr_geo_four_velocity`](@ref) its four-velocity:

```@example apex
u = kerr_geo_four_velocity(0.9, 10.0, 0.5, 0.8)
[component(1.0) for component in u]
```

## Finite-window plunge

[`kerr_geo_plunge`](@ref) builds an ``E < 1`` plunge from the constants in one of three root
classes: `"Real1"` (four real roots, three outside the horizon), `"Real2"` (four real roots,
one outside) and `"Complex"` (two real roots and a complex pair).

```@example apex
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start = :turning_point)
(plunge.OrbitClass, plunge.Status.duration.mino_time_to_horizon)
```

The orbit starts at the turning point by default. It can also start from given Mino-time
phases, `initPhases = (t0, λr0, λθ0, φ0)`, or from a position, `initial_radius` and
`initial_theta`. `time_origin = :future_horizon_v_zero` shifts ``t`` so that ``v = t + r_*``
vanishes on the future horizon.

```@example apex
(plunge.Trajectory.r(0.0), plunge.Trajectory.rstar(0.0), plunge.Trajectory.v(0.0))
```

The functions it is built from are exported as well: [`radial_roots`](@ref),
[`polar_roots`](@ref), [`classify_orbit`](@ref), [`lambda_of_r`](@ref),
[`generic_plunge_orbit`](@ref), [`generic_plunge_velocity`](@ref) and [`kerr_rstar`](@ref).

## Finite-window capture and scattering

[`kerr_geo_capture`](@ref) and [`kerr_geo_scatter`](@ref) build an ``E ≥ 1`` orbit from the
constants, `kerr_geo_capture(a, (E, Lz, Q))`, and `kerr_geo_scatter` also from the APEX
parameters of an equatorial orbit. `Formula` names the closed form of ``r(λ)``:

| `Formula` | Energy | Radial roots | ``λ = 0`` |
| :--- | :--- | :--- | :--- |
| `:hyperbolic_scatter` | ``E > 1`` | four real | closest approach |
| `:parabolic_scatter` | ``E = 1`` | three real | closest approach |
| `:hyperbolic_capture` | ``E > 1`` | two real inside the horizon and a complex pair | the future horizon |
| `:parabolic_capture` | ``E = 1`` | one real and a complex pair | the future horizon |

Constants without such an orbit give a record with `Status.supported == false`.
`polar_phase` sets the polar phase at ``λ = 0``: phase 0 is the northern turning point for
``E > 1`` and the equator, crossed northward, for ``E = 1``. For a capture, ``t`` and ``φ``
are given as increments over a finite window, since they diverge on the horizon; ``v`` and
``ψ`` are given directly.

```@example apex
scatter = kerr_geo_scatter(0.5, (1.1, 5.0, 1.0))
(scatter.Formula, scatter.Outcome, scatter.Trajectory.r(0.0))
```

[`kerr_geo_scatter_asymptotic_diagnostics`](@ref) and
[`kerr_geo_scatter_asymptotic_state`](@ref) give the asymptotic directions, the impact
parameter and the deflection of a scattered orbit, built by either interface.
