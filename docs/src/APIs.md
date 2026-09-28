# API Reference

This page documents the public interfaces intended for direct calls. Canonical input is
`(a, E, Lz, Q)`; APEX-like input `(a, p, e, x)` is converted to constants first.

## Families and members

```@docs
kerr_geodesic
KerrGeodesicFamily
kerr_geo_members
KerrGeoComponent
kerr_geo_member_class
kerr_geo_sample
```

`KerrGeodesicFamily` has one slot per class: `Stable`, `Critical` (a tuple, ordered by
role), `Plunge`, `Capture`, `Scatter`, `Trapped`; horizon- and extremal-tier members sit in
the slot of their class. `Status.case_ids` lists the cases of the constants,
`Status.member_errors` the members that could not be built, and an empty family has
`Status.supported == false` with a `Status.reason`. `kerr_geodesic` also accepts the
selection keywords of [`kerr_geo_select_component`](@ref) (`case_id`, `initial_radius`,
`radial_sign`, `endpoint_intent`), the polar keywords (`polar_sector`, `polar_phase`,
`polar_hemisphere`), `axis=:north`/`:south` for motion along the spin axis,
`reference_radius`, `initPhases` (Stable member) and `trapped_component` (`:full`,
`:outgoing`, `:incoming`). A selection keyword that matches no member, or more than one, is
an error on every tier.

The member types are aliases of `KerrGeoComponent{C}`:

| type | class | cases |
| :--- | :--- | :--- |
| `KerrGeoStableComponent` | `:stable` | A1, A2, A-H1, A-H2, A-X1, A-X2 |
| `KerrGeoCriticalComponent` | `:critical` | K1–K11 |
| `KerrGeoPlungeComponent` | `:plunge` | B1–B9, B-X1, B-X2 |
| `KerrGeoCaptureComponent` | `:capture` | C1–C12, C-X1–C-X4 |
| `KerrGeoScatterComponent` | `:scatter` | D1, D2, D-H1, D-H2, D-X1, D-X2 |
| `KerrGeoTrappedComponent` | `:trapped` | N1–N6 |

Member fields:

- `CaseId`; `Tier` (`:primary`, `:horizon`, `:extremal`); `Role` (Critical: `:on_root`,
  `:outer`, `:inner`, the side of the repeated root; otherwise `:none`); `Component`;
  `ConstantsOfMotion = (a, E, Lz, Q)`; `Roots`.
- `ReferenceZero`: the events of `λ = 0` (`lambda0_event`), of `t = φ = 0`
  (`t_phi_zero_event`, at `t_phi_zero_lambda` and radius `t_phi_zero_radius`) and of `τ = 0`
  (`tau_zero_event`, always at `λ = 0`); `lambda_regular`, where `v` and `ψ` vanish, or
  `nothing`; the polar phase (`polar_phase`, `polar_phase_convention`; `phases` for Stable
  members).
- `Domain`: `mino`, `endpoint_closed`, `endpoint_roles`, and `horizon_lambda` for members
  that end on the future horizon.
- `Trajectory`: `t, r, theta, z, phi, tau`; `rstar, v, psi` where a horizon-regular chart
  exists (`v, psi` for Trapped members with `|a| < 1`); `u, chi` for Trapped and `|a| = 1`
  members.
- `Velocity` (the Mino-time rates `ut, ur, uz, utheta, uphi, dtau_dlambda`), `Potentials`
  (`radial`, `polar_z`), `Residuals` (`radial`, `polar_z`, `normalization`).
- `Status`: `supported`, member-specific metadata (Stable: `apex`, `frequencies`,
  `precision`) and `spectral`, a `SpectralStatus`. Its `achieved` is the largest error
  estimate of the member's Chebyshev tables relative to the local size of the tabulated
  rate, `pieces` their total number, and `unresolved` is 0, because a table that misses its
  tolerance raises an error. A piece whose resolution is limited by floating-point rounding
  (of `λ` or of the rate itself) is accepted, and its estimate shows in `achieved`.

`kerr_geo_sample(m, λs)` evaluates `t, r, theta, phi, tau, ut, ur, utheta, uphi` at many
Mino times in one call. The loop runs behind a function barrier on the member's concrete
functions, which avoids the dynamic dispatch of calling `m.Trajectory.t(λ)` and the others
one by one from untyped code.

## Classes, cases and tiers

```@docs
kerr_geo_class
kerr_geo_case_class
kerr_geo_tier
kerr_geo_case_name
kerr_geo_case_symbol
kerr_geo_is_critical
kerr_geo_critical_role
kerr_geo_case_catalog
```

| name | content |
| :--- | :--- |
| `KERR_GEO_CLASSES` | the six classes: case letter, Symbol, display name, family slot, order, plot colours |
| `HORIZON_TIER_IDS` | `(:A_H1, :A_H2, :D_H1, :D_H2)`: `0 < \|a\| < 1`, `P(r₊) = 0` |
| `EXTREMAL_TIER_IDS` | `(:A_X1, :A_X2, :B_X1, :B_X2, :C_X1, :C_X2, :C_X3, :C_X4, :D_X1, :D_X2)`: `\|a\| = 1`, horizon a double or triple root |
| `CRITICAL_ROLES` | `(on_root=(:K1, :K3, :K6, :K9), outer=(:K4, :K7, :K10), inner=(:K2, :K5, :K8, :K11))` |
| `kerr_geo_case(id)` | the `KerrGeoCaseSpec` of a case: energy regime and sign, root ordering and multiplicities, allowed interval, endpoints, paired cases, formula family |

## Classification

```@docs
kerr_geo_classify
kerr_geo_components
kerr_geo_select_component
kerr_geo_classification_pipeline
kerr_geo_root_structure
```

`kerr_geo_classify` returns a `KerrGeoClassification` that keeps every admitted exterior
radial component (`KerrGeoRadialComponent`, with `CaseId`, `BroadClass`, endpoints and
`PolarSector`) before any selection; selection does not delete the other members of a
family. `kerr_geo_classification_pipeline()` gives the order of the stages: canonical
constants, metric limit, energy regime, radial degree, roots and multiplicities, allowed
components, polar admissibility, future direction, case assignment, family assembly,
optional selection.

## Member constructors

Each constructor takes `(a, E, Lz, Q)`, classifies the constants and builds one member; the
class constructors (`kerr_geo_*_component`, `kerr_geo_trapped`) also take
`(a, (E, Lz, Q))`. `kerr_geo_capture_equator_attractive(a, E, Lz)` fixes `Q = 0`, the
axis-infall constructors take `(a, E; axis)`, and `kerr_geo_trapped_case` takes the case ID
first. When the constants admit several members of the class, choose one with `case_id`.

In every member `r(λ)` is a closed form; for `|a| < 1`, `t`, `φ`, `τ`, `v`, `ψ` come from the
spectral engine.

### Stable

```@docs
kerr_geo_stable_component
```

### Critical

```@docs
kerr_geo_critical_component
kerr_geo_critical_spherical
kerr_geo_critical_homoclinic
kerr_geo_critical_plunge
```

| cases | radial motion | Mino domain | reference events |
| :--- | :--- | :--- | :--- |
| K1, K3, K6, K9 | on the repeated root `r_c` | `(-∞, ∞)` | `t`, `φ`, `τ` zero at `λ = 0` (the polar reference event); no horizon chart |
| K4 | `r_c` → apastron → `r_c` | `(-∞, ∞)` | apastron at `λ = 0`, where `t = φ = τ = 0`; `lambda_of_radius(r; branch=:outgoing/:incoming)` |
| K7, K10 | infinity → `r_c` | `(λ_∞, ∞)` | `t`, `φ`, `τ`, `v`, `ψ` zero at `λ = 0` (the reference radius) |
| K2, K5, K8, K11 | `r_c` → future horizon | `(-∞, 0]` | future horizon at `λ = 0`, where `τ = v = ψ = 0`; `t`, `φ` zero at the reference radius (`ReferenceZero.t_phi_zero_lambda < 0`); `polar_phase` is given at the reference-radius event for K2, K8, K11 and on the horizon for K5 (`ReferenceZero.polar_phase_event`) |

### Plunge

```@docs
kerr_geo_plunge_component
kerr_geo_plunge_axis_infall
```

The members of `kerr_geo_plunge_component` start at their turning point (`λ = 0`, where
`t`, `φ`, `τ` vanish) and end on the future horizon at `Domain.horizon_lambda`, where
`v = ψ = 0`. `r(λ)` comes from the shared radial models for B1, B3 and B4, from
multiplicity-specific closed forms for B2, B5 and B6, and from the elementary integrals of a
repeated root strictly inside the horizon for B7–B9. Motion along the spin axis
(`kerr_geo_plunge_axis_infall`) puts the horizon at `λ = 0` instead. The plunges from a
Critical root are K2, K5, K8, K11.

### Capture

```@docs
kerr_geo_capture_component
kerr_geo_capture_four_complex
kerr_geo_capture_vortical
kerr_geo_capture_constant_latitude
kerr_geo_capture_equator_attractive
kerr_geo_capture_axis_infall
```

| cases | radial motion | Mino domain | reference events |
| :--- | :--- | :--- | :--- |
| C1–C12 | infinity → future horizon, no turning point | `(λ_∞, 0]` | future horizon at `λ = 0`, where `τ = v = ψ = 0`; `t`, `φ` zero at `ReferenceZero.t_phi_zero_lambda` (default `λ_∞/2`, or `reference_radius`) |

C1 and C3 take `r(λ)` from the closed-form radial formulas of [`kerr_geo_capture`](@ref), C2
and C4 from their radial models, C5 (four complex radial roots, formula family FF18) from the
Jacobi model of `kerr_geo_capture_four_complex`, which accepts vortical, constant-latitude,
axis-crossing (`Lz = 0`) and exact-axis polar motion, and C6–C12 from the elementary
integrals of a repeated root strictly inside the horizon. The captures that tend to an
unstable repeated root are the Critical members K7 and K10. The domain is open at `λ_∞`
(`r → ∞`); the member's functions raise a `DomainError` there.

### Scatter

```@docs
kerr_geo_scatter_component
kerr_geo_scatter_asymptotic_diagnostics
kerr_geo_scatter_asymptotic_state
```

| cases | radial motion | Mino domain | reference events |
| :--- | :--- | :--- | :--- |
| D1 (`E = 1`), D2 (`E > 1`) | infinity → outer turning point → infinity | `(-λ_∞, λ_∞)` | turning point at `λ = 0`, where `t = φ = τ = 0`; both infinity endpoints excluded |

`Trajectory.lambda_of_radius(r; branch=:outgoing)` inverts `r(λ)` on either leg. The
asymptotic diagnostics (directions, impact data, deflection angles) are available for the
pendular and equatorial sectors. D1 shares its constants with the Plunge member B5, D2 with
B6.

### Trapped

```@docs
kerr_geo_trapped
kerr_geo_trapped_case
kerr_geo_trapped_classify
KerrGeoTrappedClassification
```

Trapped members (`E < 0`) run from the past horizon through a turning point (`λ = 0`) into
the future horizon; for `|a| = 1` they belong to the extremal family (`Tier == :extremal`).
For `0 < |a| < 1` (`kerr_geo_trapped`) the full Mino domain is `[-Λ_H, Λ_H]`
(`Domain.horizon_half_duration = Λ_H`), the outgoing half `[-Λ_H, 0]` and the incoming half
`[0, Λ_H]`. BL `t` and `φ` exclude both horizons and diverge there; the past horizon uses
the retarded pair `u`, `chi` and the future horizon the advanced pair `v`, `psi`, each zero
on its own horizon. `component=:full`, `:outgoing` or `:incoming` (`trapped_component=` in
`kerr_geodesic`) selects the Mino domain of the `Trajectory` and `Velocity` functions;
`Trajectory.full`, `Trajectory.outgoing` and `Trajectory.incoming` return the coordinates at
one `λ` as a NamedTuple. `Status.disposition_id` is the disposition label NFD0k of case Nk.

### Horizon and extremal tiers

```@docs
kerr_geo_horizon_stable
kerr_geo_horizon_scatter
kerr_geo_extremal_family
kerr_geo_extremal
KerrGeoExtremalFamily
```

`kerr_geodesic` dispatches `P(r₊) = 0` constants (`0 < |a| < 1`) to the horizon-tier
constructors (A-H1, A-H2 for `E < 1`; D-H1, D-H2 for `E ≥ 1`), and `|a| = 1` to the
extremal family. Constants are used exactly as given; APEX input `kerr_geodesic(a, p, e, x)`
takes a spin within `8 eps` of `±1` as `±1` and keeps the value passed in as
`Status.input_provenance.original`. Near-extremal spins such as `1 − 1e-8` stay
sub-extremal. At `|a| = 1` every member has `Tier == :extremal`: with
`P_H = 2E − aLz > 0` the members keep their primary case IDs, with `P_H = 0` (the horizon a
double or triple root of `R`) they are A-X1, A-X2, B-X1, B-X2, C-X1…C-X4, D-X1, D-X2. Their zero
events follow the regular chart: members that end on the future horizon (B, C, K2, K5, K8,
K11) have `λ = 0` and `v = ψ = 0` there, with `t`, `φ` fixed by that chart; A1 has
`λ = 0` and `t = φ = 0` at the inner turning point; members that reach no horizon have
`lambda_regular = nothing`.
`kerr_geo_extremal_family` returns the `KerrGeoExtremalFamily` (`MetricLimit`,
`ConstantsOfMotion`, `Classification`, `Members`, `Status`) directly.

The signed-spin reflection `(a, Lz, φ) → (−a, −Lz, −φ)` at fixed `E`, `Q`, Mino origin and
polar phase maps `a = 1` members onto `a = −1` members: `r`, `z`, `t`, `τ`, `rstar`, `u`,
`v` are even and `phi`, `psi`, `chi` odd.

## APEX and finite-window APIs

These functions work with APEX parameters `(a, p, e, x)` or build one finite-window orbit
from the constants, and return their own record types rather than `KerrGeoComponent`
members.

```@docs
kerr_geo_stable
KerrGeoStable
kerr_geo_orbit_type_metadata
kerr_geo_orbit_type
kerr_geo_four_velocity
```

| function | purpose |
| :--- | :--- |
| `kerr_geo_constants_of_motion(a,p,e,x)` | `Dict("E"=>…, "Lz"=>…, "Q"=>…)` |
| `kerr_geo_frequencies(a,p,e,x; Time="Mino")` | Mino, `"BoyerLindquist"` or `"Proper"` frequencies |
| `kerr_geo_orbit(a,p,e,x; initPhases=…)` | dictionary-style orbit with the key `"Stability"`, for Stable orbits, constant-radius Critical orbits with `E < 1` (ISCO/ISSO, unstable circular and spherical orbits: `ϒr = 0` and, when unstable, `"RadialLyapunovExponent"` = √(R''/2) in Mino time); any other input is an `ArgumentError` pointing to `kerr_geodesic(a, p, e, x)` |
| `kerr_geo_radial_roots(a,p,e,x)`, `kerr_geo_polar_roots(a,p,e,x)` | radial and polar roots of an APEX point |

`kerr_geo_orbit_type_metadata(a, p, e, x)` gives `family` ("Stable", "Critical",
"Plunge", "Capture", "Scatter"; unstable circular orbits, the ISCO/ISSO and orbits on the
separatrix are "Critical"), `outcome` (the class Symbol), `labels = [family, shape,
inclination]`, `stability` ("Stable", "MarginallyStable", "Unstable", "NotApplicable") and
`energy_regime` ("Elliptic", "Parabolic", "Hyperbolic"). Inputs with
`|p − separatrix_p| ≤ 1e-15` are evaluated at `separatrix_p` (`at_separatrix`), and orbits
within `1e-12` of the separatrix are "Critical". A constant-radius orbit with
`E ≥ 1` (K6, K9) is `kerr_geodesic(a, p, e, x).Critical[1]`.

### Finite-window plunge

```@docs
kerr_geo_plunge
KerrGeoPlunge
kerr_rstar
```

| function | purpose |
| :--- | :--- |
| `radial_roots(a,E,L,Q)` | radial-potential roots |
| `polar_roots(a,E,L,Q)` | polar roots |
| `classify_orbit(a,E,L,Q; atol=1e-15)` | `(roots, class)` with class `"Complex"`, `"Real1"` or `"Real2"`; any other root structure raises an error |
| `lambda_of_r(a,E,L,Q)` | `(Λ_end, Λ_H, λ_of_r)`: the Mino times from the turning point to the inner end of the radial range and to the outer horizon, and the map `r ↦ λ` |
| `generic_plunge_orbit(a,E,L,Q; initPhases=…)` | low-level plunge trajectory callables |
| `generic_plunge_velocity(a,E,L,Q; initPhase=…)` | low-level plunge velocity callables |

`KerrGeoPlunge.OrbitClass` carries the root class (`Complex`, `Real1`, `Real2`);
`Status.duration` records the finite Mino time to the horizon.

### Finite-window capture and scatter

```@docs
kerr_geo_capture
KerrGeoCapture
kerr_geo_scatter
KerrGeoScatter
```

`kerr_geo_scatter` builds the `E ≥ 1` finite-window scatter orbit from constants
(`input=:constants`, or `kerr_geo_scatter(a, (E, Lz, Q))`) or from the APEX parameters of an
equatorial orbit. `kerr_geo_capture` builds the capture orbit from constants
(`input=:constants`, or `kerr_geo_capture(a, (E, Lz, Q))`); a capture orbit has no
periapsis, so it has no APEX parameters. The result's `Formula` names the closed form used:

| `Formula` | energy | radial roots | reference zero |
| :--- | :--- | :--- | :--- |
| `:hyperbolic_scatter` | E > 1 | four real | closest approach |
| `:parabolic_scatter` | E = 1 | three real | closest approach |
| `:hyperbolic_capture` | E > 1 | two real inside the horizon + complex pair | future horizon (regular v, ψ = 0) |
| `:parabolic_capture` | E = 1 | one real + complex pair | future horizon (regular v, ψ = 0) |

`Outcome` is the outcome at infinity (`:scatter`, `:capture`, …) and `EnergyRegime` the
energy regime (`:parabolic`, `:hyperbolic`, …). Constants without such a component return
a record with `Status.supported == false`.
`polar_phase` (default 0) fixes the polar phase at the reference event (closest approach for
scatter, the horizon for capture); phase 0 is the northern polar turning point for `E > 1`
and the equator, crossed northward, for `E = 1`. Boyer–Lindquist t and φ of capture
orbits diverge on the horizon, so they are given as finite-window increments, while the
ingoing coordinates v and ψ are given directly.

## Orbital landmarks

```@docs
kerr_geo_separatrix
kerr_geo_isco
kerr_geo_ibso
kerr_geo_isso
```

## Metric functions

| function | purpose |
| :--- | :--- |
| `kerr_horizons(a)` | `(rplus, rminus)` |
| `kerr_delta(a, r)` | `Δ = r² − 2r + a²` |
| `kerr_metric_limit(a)` | `:schwarzschild`, `:subextremal`, `:near_extremal`, `:extremal` |
| `kerr_energy_regime(E)` | `:elliptic`, `:parabolic`, `:hyperbolic`: the sign of `E² − 1` |
| `kerr_energy_sign(E)` | `-1` for `E < 0`, else `+1` |
| `kerr_radial_momentum(a,E,Lz,r)` | `P(r) = E(r² + a²) − a Lz` |
| `kerr_radial_potential(a,E,Lz,Q,r)` | `R(r)` |
| `kerr_radial_polynomial`, `kerr_radial_coefficients`, `kerr_radial_derivatives` | `R` as a polynomial, its coefficients and derivatives |
| `kerr_polar_z_potential(a,E,Lz,Q,z)` | `Θ(z) = (dz/dλ)²`, `z = cos θ` |
| `kerr_polar_theta_potential(a,E,Lz,Q,θ)` | `Θ(θ)` |
| `kerr_polar_admissibility`, `kerr_polar_sector_candidates` | polar admissibility and the possible polar sectors of the constants |
| `kerr_axis_carter_q(a,E)` | `Q = a²(1 − E²)` of motion along the spin axis |

## Diagnostics

```@docs
kerr_geo_diagnose
```
