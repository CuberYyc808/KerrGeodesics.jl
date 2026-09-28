# Examples

Each section builds one kind of orbit and reads off what is particular to it. The notebook
`example/KerrGeodesics_56_Orbit_Catalog.ipynb` holds one set of constants for each of the
56 cases, with plots and the animated gallery shown on the [home page](index.md).

```@setup ex
using KerrGeodesics
```

## A stable orbit

```@example ex
kg = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
stable = kg.Stable
(stable.CaseId, stable.Trajectory.r(0.0), stable.Trajectory.z(0.0))
```

With zero initial phases, ``λ = 0`` is at periapsis, ``r = p/(1 + e)``, and at the northern
polar turning point. `initPhases` moves the start along the orbit:

```@example ex
shifted = kerr_geodesic(0.9, 10.0, 0.5, 0.8; initPhases = (0.0, π, 0.0, 0.0)).Stable
shifted.Trajectory.r(0.0)                       # apoapsis, p/(1 - e)
```

A circular equatorial orbit has ``e = 0`` and ``x = ±1``:

```@example ex
circular = kerr_geodesic(0.9, 8.0, 0.0, 1.0).Stable
(circular.CaseId, circular.Trajectory.r(3.0), circular.Trajectory.theta(3.0))
```

## Critical orbits

These constants have an unstable spherical orbit at ``a = 0.7``. The family holds its three
Critical members, ordered by role:

```@example ex
kg = kerr_geodesic(0.7, (0.9171300256198305, 2.2591913519439517, 2.898994491984013))
[(member.CaseId, member.Role, member.Status.name) for member in kg.Critical]
```

The homoclinic orbit K4 leaves the spherical orbit, turns at apoapsis at ``λ = 0`` and comes
back, so every radius off the root is reached twice. `lambda_of_radius` takes the branch:

```@example ex
homoclinic = kg.Critical[2]
r1 = homoclinic.Trajectory.r(1.0)
(homoclinic.Trajectory.r(0.0),
 homoclinic.Trajectory.lambda_of_radius(r1; branch = :incoming),
 homoclinic.Trajectory.lambda_of_radius(r1; branch = :outgoing))
```

The whirling orbit K5 leaves the spherical orbit inward and crosses the horizon at
``λ = 0``, where ``v`` vanishes:

```@example ex
whirl = kg.Critical[3]
(whirl.Domain.mino, whirl.Trajectory.v(0.0))
```

The ISCO or ISSO is the Critical member K1, at ``p =`` [`kerr_geo_isso`](@ref)`(a, x)` and
``e = 0``:

```@example ex
isso = kerr_geo_isso(0.9, 0.8)
kerr_geodesic(0.9, isso, 0.0, 0.8).Critical[1].CaseId
```

A single Critical member can also be built on its own:

```@example ex
from_infinity = kerr_geo_critical_component(0.7, 1.0, -0.7, 16.0; case_id = :K7)
from_infinity.Domain.endpoint_roles
```

## A plunge

```@example ex
plunge = kerr_geodesic(0.9, (0.94, 0.1, 12.0)).Plunge
λH = plunge.Domain.horizon_lambda
(plunge.CaseId, plunge.Trajectory.r(0.0), plunge.Trajectory.r(λH))
```

The plunge starts at its turning point and crosses the horizon at ``λ_H``, where ``v`` and
``ψ`` vanish and Boyer–Lindquist ``t`` diverges:

```@example ex
(plunge.Trajectory.v(λH), plunge.Trajectory.psi(λH), plunge.Trajectory.t(0.999λH))
```

## A capture

```@example ex
capture = kerr_geo_capture_component(0.9, 1.1, 0.5, 3.0)
(capture.CaseId, capture.Domain.mino)
```

The capture comes in from infinity at ``λ_∞ < 0`` and crosses the horizon at ``λ = 0``. ``t``
and ``φ`` vanish halfway in Mino time, at `ReferenceZero.t_phi_zero_lambda`:

```@example ex
(capture.ReferenceZero.t_phi_zero_lambda, capture.ReferenceZero.t_phi_zero_radius,
 capture.Trajectory.v(0.0))
```

With four complex radial roots and ``Q < 0`` the orbit is C5, and its polar motion is
vortical: it stays in one hemisphere.

```@example ex
vortical = kerr_geo_capture_four_complex(0.9, 1.8, 0.2, -1.0; polar_hemisphere = :north)
λmid = vortical.Domain.mino[1] / 2
(vortical.CaseId, vortical.Status.polar.sector, vortical.Trajectory.z(λmid))
```

With ``L_z = 0`` and ``a^2(1 - E^2) < Q < 0``, the same constants can also describe motion
over the poles, and `polar_sector` chooses it:

```@example ex
a, E = 0.9, 1.5
over_the_poles = kerr_geo_capture_four_complex(a, E, 0.0, a^2 * (1 - E^2) + 1e-5;
                                               polar_sector = :axis_crossing)
over_the_poles.Status.polar.sector
```

Motion along the spin axis takes the axis explicitly:

```@example ex
along_axis = kerr_geo_capture_axis_infall(0.5, 1.2; axis = :north)
(along_axis.CaseId, along_axis.Trajectory.theta(along_axis.Domain.mino[1] / 2))
```

## A scattered orbit

```@example ex
kg = kerr_geodesic(0.5, (1.1, 5.0, 1.0))
scatter = kg.Scatter
(kg.Status.case_ids, scatter.Domain.mino, scatter.ReferenceZero.t_phi_zero_radius)
```

The orbit turns at closest approach at ``λ = 0``; each radius is reached once on the way in
and once on the way out. The asymptotic data give the directions of approach and escape and
the deflection:

```@example ex
asymptotics = kerr_geo_scatter_asymptotic_diagnostics(scatter)
(asymptotics.azimuthal_deflection.value, asymptotics.deflection_angle_3d.value,
 asymptotics.impact_magnitude.value)
```

## A trapped orbit

With ``E < 0`` the orbit lives inside the ergoregion. It leaves the past horizon at
``-λ_H``, turns at ``λ = 0`` and crosses the future horizon at ``λ_H``:

```@example ex
trapped = kerr_geo_trapped(0.9, -0.8, -4.0, 1.0)
ΛH = trapped.Domain.horizon_half_duration
(trapped.CaseId, trapped.Domain.mino, trapped.Trajectory.full(0.0).r)
```

``t`` and ``φ`` diverge on both horizons. The outgoing coordinates ``u``, ``χ`` vanish on the
past horizon and the ingoing ``v``, ``ψ`` on the future one:

```@example ex
(trapped.Trajectory.u(-ΛH), trapped.Trajectory.v(ΛH))
```

`component = :incoming` keeps the half after the turning point; in `kerr_geodesic` the
keyword is `trapped_component`:

```@example ex
incoming = kerr_geodesic(0.9, (-0.8, -4.0, 1.0); trapped_component = :incoming).Trapped
incoming.Domain.mino
```

## Horizon and extremal members

Constants with ``P(r_+) = 0`` make the horizon a root of ``R``:

```@example ex
a = 0.9
rplus = 1 + sqrt(1 - a^2)
island = kerr_geodesic(a, (0.94, 2rplus * 0.94 / a, 0.0)).Stable
(island.CaseId, island.Tier, kerr_geo_case_name(island.CaseId))
```

At ``a = ±1`` every member belongs to the extremal tier. Reflecting the spin and ``L_z``
leaves ``r`` unchanged and reverses ``φ``:

```@example ex
plus = kerr_geodesic(1.0, (1.2, 2.0, 14.0)).Scatter
minus = kerr_geodesic(-1.0, (1.2, -2.0, 14.0)).Scatter
λ = 0.2 * plus.Domain.mino[2]
(plus.Tier, plus.Trajectory.r(λ) - minus.Trajectory.r(λ),
 plus.Trajectory.phi(λ) + minus.Trajectory.phi(λ))
```

With ``P(r_+) = 2E - aL_z = 0`` at ``a = 1`` the horizon is a double or triple root, and the
members have their own IDs:

```@example ex
kerr_geodesic(1.0, (0.8, 1.6, 1.0)).Stable.CaseId
```
