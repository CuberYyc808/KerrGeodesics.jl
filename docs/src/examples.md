# Examples

## The 56-orbit catalogue

The notebook `example/KerrGeodesics_56_Orbit_Catalog.ipynb` (with its helpers in
`example/kerr_geodesics_56_case_support.jl`) holds one fixed set of constants for every
catalogue orbit, keyed by its case ID (A1–A2, K1–K11, B1–B9, C1–C12, D1–D2, N1–N6) or, for
the tier members, by their tier name (A-H1, A-H2, A-X1, A-X2, B-X1, B-X2, C-X1…C-X4, D-H1,
D-H2, D-X1, D-X2). The 56 orbits split into Stable 6, Critical 11, Plunge 11, Capture 16,
Scatter 6 and Trapped 6. For each one it calls `kerr_geodesic`, shows the classification,
picks the member from its class slot and plots the trajectory; its last cell animates all
56 orbits in the style of `example/animations/showcase_all.gif`, the grid of all 56 shown in
the README.

## Stable orbit

```julia
using KerrGeodesics

family = kerr_geodesic(0.9, 10.0, 0.5, 0.8)     # APEX-like (a, p, e, x)
family.Status.case_ids                          # (:A1, :B1)

m = family.Stable                               # KerrGeoStableComponent, case A1
t = m.Trajectory.t
r = m.Trajectory.r
theta = m.Trajectory.theta
phi = m.Trajectory.phi
m.Status.apex                                   # (a, p, e, x) from the turning points
m.Status.frequencies                            # Mino-time frequencies
```

`λ = 0` is the periapsis and the northern polar turning point; the keyword
`initPhases = (qt0, qr0, qθ0, qφ0)` shifts the phases as in `kerr_geo_orbit`. The same
member from constants is `kerr_geo_stable_component(a, E, Lz, Q)`.

The reference implementation `kerr_geo_stable` builds the same orbit from APEX input:

```julia
stable = kerr_geo_stable(0.9, 10.0, 0.5, 0.8; initPhases=(0.0, 0.0, 0.0, 0.0))
stable.Trajectory.r(0.0)                        # p/(1+e): periapsis
stable.Frequencies
```

Circular equatorial orbits have constant `r = p` and `θ = π/2`. The class of an APEX point
is

```julia
metadata = kerr_geo_orbit_type_metadata(0.9, 10.0, 0.5, 0.8)
metadata.family                                 # "Stable"
metadata.labels                                 # ["Stable", "Eccentric", "Inclined"]
metadata.stability                              # "Stable"
metadata.at_separatrix                          # false
```

## Critical orbits

Constants with an unstable repeated root admit one Critical member per role. At `a = 0.7`:

```julia
family = kerr_geodesic(0.7, (0.9171300256198305, 2.2591913519439517, 2.898994491984013))
Tuple(m.CaseId for m in family.Critical)        # (:K3, :K4, :K5)
Tuple(m.Role for m in family.Critical)          # (:on_root, :outer, :inner)
Tuple(m.Status.name for m in family.Critical)   # (:unstable_spherical, :homoclinic, :whirling)

homoclinic = family.Critical[2]                 # K4: r_c -> apastron (λ = 0) -> r_c
homoclinic.Domain.mino                          # (-Inf, Inf)
homoclinic.Trajectory.r(0.0)                    # the apastron
homoclinic.Trajectory.lambda_of_radius(homoclinic.Trajectory.r(1.0); branch=:incoming)  # 1.0

whirling = family.Critical[3]                   # K5: r_c -> future horizon at λ = 0
whirling.Trajectory.v(0.0)                      # 0.0
```

Single members come from `kerr_geo_critical_spherical`, `kerr_geo_critical_homoclinic`,
`kerr_geo_critical_plunge`, or `kerr_geo_critical_component(a, E, Lz, Q; case_id)`:

```julia
k7 = kerr_geo_critical_component(0.7, 1.0, -0.7, 16.0; case_id=:K7)   # E = 1, from infinity
k7.Domain.endpoint_roles                        # (:past_infinity, :future_repeated_root_asymptote)
```

The ISCO/ISSO (`p = kerr_geo_isso(a, x)`, `e = 0`) is the Critical member K1:

```julia
isso = kerr_geo_isso(0.9, 0.8)
kerr_geodesic(0.9, isso, 0.0, 0.8).Critical[1].CaseId     # :K1
```

## Plunge from constants

```julia
family = kerr_geodesic(0.9, (0.94, 0.1, 12.0))
m = family.Plunge                               # B4: two real roots and a complex pair
m.Domain.mino                                   # (0.0, λ_H): turning point to future horizon
λH = m.Domain.horizon_lambda
m.Trajectory.r(λH)                              # r₊
m.Trajectory.v(λH), m.Trajectory.psi(λH)        # (0.0, 0.0)
m.Trajectory.t(0.5λH)                           # BL t, defined on [0, λ_H)
family.RootClass
```

### Finite-window plunge API

The reference plunge API takes the same constants:

```julia
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0; radial_start=:turning_point)

plunge.Trajectory.r(0.0)
plunge.Trajectory.u(0.0)
plunge.Trajectory.v(0.0)
plunge.Trajectory.rstar(0.0)
plunge.Status.duration.mino_time_to_horizon     # finite Mino time to the horizon
plunge.OrbitClass                               # "Complex"
```

It can be initialized with explicit Mino-time phases,

```julia
plunge = kerr_geo_plunge(a, E, Lz, Q;
                         initPhases=(t0, lambda_r0, lambda_theta0, phi0))
```

or with initial positions:

```julia
plunge = kerr_geo_plunge(a, E, Lz, Q;
                         initial_radius=r0,
                         initial_theta=theta0)
```

`initial_theta` defaults to `π/2`, the equator. With
`time_origin=:future_horizon_v_zero` the additive time origin is shifted so that the
advanced time on the future horizon is $v_H=0$:

```julia
plunge = kerr_geo_plunge(a, E, Lz, Q; time_origin=:future_horizon_v_zero)
plunge.Status.horizon_time_shift
plunge.Status.horizon_v_anchor_method
```

This does not make the retarded time finite at the future horizon; $u=v-2r_*$ diverges as
$r_*\to-\infty$. For the Real2 root class (four real roots, one of them outside `r₊`),
`plunge.Trajectory.v_rstar_series(rstar)` and `u_rstar_series(rstar)` evaluate the same
convention near the horizon as functions of `rstar`; for the other root classes they return
`NaN`.

## Capture

```julia
c3 = kerr_geo_capture_component(0.9, 1.1, 0.5, 3.0)      # C3
c3.Domain.mino                                  # (λ_∞, 0.0): infinity to future horizon
λ = 0.5 * c3.Domain.mino[1]
c3.Trajectory.r(λ), c3.Trajectory.theta(λ)
c3.Trajectory.v(0.0)                            # 0.0 on the horizon
c3.ReferenceZero.t_phi_zero_lambda              # where t, φ vanish (τ vanishes at λ = 0)
```

Four complex radial roots (C5, formula family FF18):

```julia
c5 = kerr_geo_capture_four_complex(0.9, 1.8, 0.2, -1.0;
                                   polar_hemisphere=:north, polar_phase=0.0)
c5.CaseId                                       # :C5
c5.Status.formula_family                        # :FF18
lambda_mid = 0.5 * c5.Domain.mino[1]
c5.Trajectory.r(lambda_mid)
c5.Trajectory.theta(lambda_mid)
```

For constant-latitude C5 the constants select the double polar root; exact axis constants
use `Lz = 0` and `Q = a²(1 − E²)`, and `polar_hemisphere=:north` or `:south` fixes the
branch. The negative-`Q` axis-crossing branch (`Lz = 0`, `a²(1 − E²) < Q < 0`) needs the
polar sector explicitly:

```julia
a = 0.9
E = 1.5
axis_crossing = kerr_geo_capture_four_complex(
    a, E, 0.0, a^2 * (1 - E^2) + 1e-5;
    polar_sector=:axis_crossing,
    polar_hemisphere=:north)

axis_crossing.CaseId                     # :C5
axis_crossing.Status.polar.sector        # :axis_crossing
```

At `Q = 0` the same polar sector uses the elementary `sech` separatrix and belongs to C3:

```julia
q_zero_crossing = kerr_geo_capture_component(
    0.9, 1.5, 0.0, 0.0;
    polar_sector=:axis_crossing)

q_zero_crossing.CaseId                    # :C3
q_zero_crossing.Status.polar.formula_kind # :hyperbolic_sech_axis_crossing
```

Capture along the spin axis needs an explicit axis:

```julia
axis_capture = kerr_geo_capture_axis_infall(0.5, 1.2; axis=:north)
axis_capture.CaseId                      # :C3
axis_capture.Domain.mino
axis_capture.Trajectory.v(0.0)           # zero at the future horizon
```

## Scatter

```julia
family = kerr_geodesic(0.5, (1.1, 5.0, 1.0))
family.Status.case_ids                   # (:B6, :D2)
d2 = family.Scatter
d2.Domain.mino                           # (-λ_∞, λ_∞); turning point at λ = 0
d2.ReferenceZero.t_phi_zero_radius      # the turning radius
d2.Trajectory.lambda_of_radius(2 * d2.ReferenceZero.t_phi_zero_radius; branch=:outgoing)   # > 0
kerr_geo_scatter_asymptotic_diagnostics(d2)
```

## Trapped orbit

The full Class N trajectory (these constants are case N4) has its radial turning event at
`λ = 0`:

```julia
trapped = kerr_geo_trapped(0.9, -0.8, -4.0, 1.0)
trapped.CaseId                                # :N4
LambdaH = trapped.Domain.horizon_half_duration

turn = trapped.Trajectory.full(0.0)           # (t, r, θ, φ, z, τ) at the turning point
past_u = trapped.Trajectory.u(-LambdaH)       # zero at the past horizon
future_v = trapped.Trajectory.v(LambdaH)      # zero at the future horizon
```

BL `t` and `phi` are evaluated only for `abs(λ) < LambdaH`. One half is selected with
`component` (or `trapped_component` in `kerr_geodesic`):

```julia
incoming = kerr_geo_trapped(0.9, -0.8, -4.0, 1.0; component=:incoming)
incoming.Domain.mino                          # (0.0, LambdaH)
incoming.Trajectory.r(0.4 * LambdaH)

family = kerr_geodesic(0.9, (-0.8, -4.0, 1.0); trapped_component=:incoming)
family.Trapped.CaseId                         # :N4
family.Trapped.Status.disposition_id          # :NFD04
```

## Horizon and extremal tiers

Constants with `P(r₊) = 0` have their own members:

```julia
a = 0.9
rplus = 1 + sqrt(1 - a^2)
island = kerr_geodesic(a, (0.94, 2rplus * 0.94 / a, 0.0)).Stable
island.CaseId, island.Tier                    # (:A_H1, :horizon)
kerr_geo_case_name(island.CaseId)             # "A-H1"
```

Exact `|a| = 1` uses its own formulas; every member has `Tier == :extremal`. These
mirrored constants give the same radial and polar motion with opposite azimuth:

```julia
plus = kerr_geodesic(1.0, (1.2, 2.0, 14.0))
minus = kerr_geodesic(-1.0, (1.2, -2.0, 14.0))

plus.Status.metric_limit       # :extremal_plus
minus.Status.metric_limit      # :extremal_minus
plus.Status.case_ids           # (:B6, :D2)
minus.Status.case_ids          # (:B6, :D2)

plus_scatter = plus.Scatter
minus_scatter = minus.Scatter
lambda_probe = 0.2 * plus_scatter.Domain.mino[2]

plus_scatter.Trajectory.r(lambda_probe)       # equal
minus_scatter.Trajectory.r(lambda_probe)
plus_scatter.Trajectory.phi(lambda_probe)     # opposite
minus_scatter.Trajectory.phi(lambda_probe)
```

With the horizon a root of `R` at `a = 1` the members have their X names, e.g.
`kerr_geodesic(1.0, (0.8, 1.6, 1.0)).Stable.CaseId == :A_X1`.
`kerr_geo_extremal_family` returns the exact-extremal family object directly. A
near-extremal value such as `1 - 1e-8` is not treated as exact.
