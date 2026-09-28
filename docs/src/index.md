# KerrGeodesics.jl

KerrGeodesics.jl builds timelike geodesics of the Kerr spacetime in units with
$G=c=M=1$, Boyer–Lindquist coordinates and Mino time $\lambda$
($d\tau/d\lambda=\Sigma$). Given the spin `a` and the constants of motion
`(E, Lz, Q)`, it classifies the radial motion outside the outer horizon and builds a member
for each radial range the constants allow, as functions of $\lambda$:
`t, r, θ (z = cos θ), φ, τ`, and the horizon-regular coordinates `v, ψ` wherever the orbit
reaches a horizon.

## Installation

```julia
using Pkg
Pkg.add("KerrGeodesics")
```

## One constructor, six classes

The unified constructor is [`kerr_geodesic`](@ref). It returns a
[`KerrGeodesicFamily`](@ref) with one slot per class:

| Class | Symbol | Cases | Family slot | Motion |
| :--- | :--- | :--- | :--- | :--- |
| Stable | `:stable` | A1, A2 | `Stable` | libration between two turning points; stable circular or spherical orbit |
| Critical | `:critical` | K1–K11 | `Critical` (tuple) | on, or asymptotic to, an unstable or marginal repeated root `r_c > r₊` |
| Plunge | `:plunge` | B1–B9 | `Plunge` | from a turning point into the future horizon |
| Capture | `:capture` | C1–C12 | `Capture` | from infinity into the future horizon (`E ≥ 1`) |
| Scatter | `:scatter` | D1, D2 | `Scatter` | from infinity through a turning point back to infinity |
| Trapped | `:trapped` | N1–N6 | `Trapped` | `E < 0`: past horizon → turning point → future horizon, inside the ergoregion |

The Critical members of one repeated root are ordered by their role: on the root
(`:on_root`: K1 ISCO/ISSO, K3/K6/K9 unstable circular or spherical orbits), on its outer
side (`:outer`: K4 homoclinic, K7/K10 from infinity) and on its inner side into the
future horizon (`:inner`: K2, K5, K8, K11).

Two degenerate situations have members outside the primary numbering, in the slot of their
class: the **horizon tier** (`0 < |a| < 1` with `P(r₊) = 0`, the horizon a root of `R`:
A-H1, A-H2, D-H1, D-H2) and the **extremal tier** (`|a| = 1` with the horizon a double or
triple root of `R`: A-X1, A-X2, B-X1, B-X2, C-X1…C-X4, D-X1, D-X2). A member's tier is
`m.Tier` (`:primary`, `:horizon` or `:extremal`); every member at `|a| = 1` has
`m.Tier == :extremal`, including those with primary IDs, while `kerr_geo_tier(id)` gives the
tier of an ID. In code the IDs are written `:A_H1`, …, and
`kerr_geo_case_name(:A_H1) == "A-H1"`.

```julia
using KerrGeodesics

family = kerr_geodesic(0.9, (0.94, 0.1, 12.0)) # constants of motion: a, (E, Lz, Q)
family_apex = kerr_geodesic(0.9, 10.0, 0.5, 0.8)   # APEX parameters: a, p, e, x

family.Status.case_ids          # the cases of these constants
kerr_geo_members(family)        # every member, in class order
```

## Members

Every member is a [`KerrGeoComponent`](@ref)`{C}`, `C` being its class
(`KerrGeoStableComponent`, `KerrGeoCriticalComponent`, `KerrGeoPlungeComponent`,
`KerrGeoCaptureComponent`, `KerrGeoScatterComponent`, `KerrGeoTrappedComponent`), with the
same fields for all classes: `CaseId`, `Tier`, `Role`, `Component`, `ConstantsOfMotion`,
`Roots`, `ReferenceZero`, `Domain`, `Trajectory`, `Velocity`, `Potentials`, `Residuals`,
`Status`. The fields of `m.Trajectory` and `m.Velocity` are functions of `λ`:

```julia
m = family.Plunge
m.CaseId                     # :B4
m.Domain.mino                # Mino-time domain: (0, λ_H), turning point to future horizon
m.Trajectory.r(0.1)
m.Trajectory.theta(0.1)
m.Trajectory.v(m.Domain.horizon_lambda)    # 0: v and ψ vanish on the future horizon
m.Velocity.ut(0.1) / m.Velocity.dtau_dlambda(0.1)   # dt/dτ
kerr_geo_member_class(m)     # :plunge
```

## Parameter conventions

The main input is the spin and the constants of motion, `(a, E, Lz, Q)`; an orbit with a
periapsis can also be given by its APEX parameters `(a, p, e, x)`. Both are four numbers, so
the constants go in a tuple (or NamedTuple): `kerr_geodesic(a, (E, Lz, Q))` reads constants,
`kerr_geodesic(a, p, e, x)` reads APEX parameters, converts them to constants and runs the
same classifier. Trapped members (`E < 0`) have no APEX
parametrization and are built from constants; exact `|a| = 1` works with either input:

```julia
kerr_geodesic(0.9, (-0.8, -4.0, 1.0)).Trapped      # N4
kerr_geodesic(1.0, (1.2, 2.0, 14.0))                # exact a = 1: B6 and D2
```

The polar motion is classified separately (`pendular`, `vortical`, `equatorial`,
`equator_attractive`, `constant_latitude`, `axis_crossing`, `axis_constant`). Where the
constants leave a choice open, keywords select it: `polar_sector=` where more than one sector
is possible, `polar_hemisphere=` (default `:north`) for motion confined to one hemisphere,
and `polar_phase=` (default `0`) for the polar phase at the reference event; what phase 0
means in each sector is recorded in `m.ReferenceZero.polar_phase_convention`.

## Time coordinates

Boyer–Lindquist `t` and `φ` diverge on a horizon. Members that reach a horizon also give
`rstar` and the ingoing coordinates $v=t+r_*$ and `psi`, which are finite there and vanish
at `m.ReferenceZero.lambda_regular`; Trapped members (and the exact `|a| = 1` members) also
give the outgoing pair `u`, `chi` for the past horizon. `τ` vanishes at `λ = 0`, and `t`,
`φ` at `m.ReferenceZero.t_phi_zero_lambda`; members that end on the future horizon from a
repeated root or from infinity (K2, K5, K8, K11, C) put the horizon at `λ = 0`.

`t, φ, τ, v, ψ` are Mino-time integrals of rates that split into radial and polar parts;
each part is stored as adaptive piecewise Chebyshev series, with the logarithmic horizon
terms of `r*` and `φ_H` added in closed form. The polar series, and the radial series of
Stable members (which give the frequencies), are fitted when the member is built; the other
radial series on the first evaluation of `t, φ, τ, v` or `ψ`. At `|a| = 1` the radial parts
are instead assembled from the closed-form Mino-time integrals of the radial models, with a
regular series at a simple horizon. `r(λ)` and `z(λ)` are closed forms (Jacobi elliptic or
elementary functions).

## APEX and finite-window APIs

The APEX API (`kerr_geo_orbit`, [`kerr_geo_stable`](@ref), `kerr_geo_frequencies`,
[`kerr_geo_orbit_type_metadata`](@ref)) and the finite-window APIs
([`kerr_geo_plunge`](@ref), [`kerr_geo_capture`](@ref), [`kerr_geo_scatter`](@ref)) are
independent reference implementations of the same geodesics. `kerr_geo_orbit` builds Stable
orbits, constant-radius Critical orbits with `E < 1` (ISCO/ISSO, unstable circular and
spherical orbits), and throws an `ArgumentError` for anything else, which
`kerr_geodesic(a, p, e, x)` covers; the class of an APEX point is
`kerr_geo_orbit_type_metadata(a, p, e, x).family`.

## Accuracy limits

These limits come from Float64 itself, not from the method. `m.Status.precision` gives the
numbers for Stable members, and `m.Status.spectral` the accuracy reached by the Chebyshev
tables of every member.

- **Polar phase.** `z(λ)` is a Jacobi function of `u = u₀ + ωλ`; rounding `u` leaves an
  absolute phase error of about `ε|ωλ|` (`ε` the machine epsilon), which reaches order one
  at `|λ| ~ 1/(ω ε)`.
- **Nearly parabolic stable orbits.** One radial period advances `t` by
  `m.Status.precision.t_radial_period`; after it `t` is known to
  `m.Status.precision.t_ulp_per_period` at best (about `10³ M` for `|E − 1| ~ 10⁻¹³`).
- **Near the spin axis.** On the narrow azimuthal spike of an orbit with small `Lz`,
  `φ(λ)` carries the error `|dφ/dλ| · ulp(λ)` from the rounding of `λ` (about `10⁻³` rad
  for `Lz = 10⁻¹²`); `Δφ` over a period from any point off the spike is accurate to
  `~10⁻¹⁵`.
- **Far from the hole.** For `E > 1`, `r(λ) ≈ 1/(√(E² − 1) |λ − λ_∞|)`, so `λ` resolves
  radii only up to `~1/(√(E² − 1) ulp(λ_∞))` (`~10¹⁶`); farther out use the radius-based
  increments `m.Trajectory.radial_*_increment(r1, r2)`, which take radii directly.

## Source layout

`src/core` (Chebyshev tools, metric functions, the polar and radial coordinate engines) →
`src/classify` (class table, case table, classifiers) → `src/models` (radial and polar
motion models) → `src/members` (the member type and its shared assembly; one directory per
class, plus the exact-extremal members; `kerr_geo_stable` and the finite-window capture and
scatter APIs sit with their class) → `src/interfaces` (APEX constants, frequencies and
`kerr_geo_orbit`; the finite-window plunge API) → `src/family` (`kerr_geodesic`) →
`src/Diagnostics.jl` (`kerr_geo_diagnose`).
