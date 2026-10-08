# Orbit classes

## What decides the motion

In Mino time ``λ``, defined by ``dτ/dλ = Σ = r^2 + a^2\cos^2θ``, the radial and polar
motions of a Kerr geodesic separate:

```math
\left(\frac{dr}{dλ}\right)^2 = R(r), \qquad \left(\frac{dz}{dλ}\right)^2 = Θ(z), \qquad z = \cos θ,
```

```math
R(r) = \bigl[E(r^2+a^2) - aL_z\bigr]^2 - Δ\,\bigl[r^2 + (L_z - aE)^2 + Q\bigr],
\qquad Δ = r^2 - 2r + a^2,
```

```math
Θ(z) = Q\,(1 - z^2) - z^2\bigl[L_z^2 + a^2(1 - E^2)(1 - z^2)\bigr].
```

An orbit can only be where ``R ≥ 0``. The leading term of ``R`` is ``(E^2 - 1)\,r^4``, and
its sign sets the energy regime, returned by [`kerr_energy_regime`](@ref):

| Regime | Energy | At large ``r`` |
| :--- | :--- | :--- |
| `:elliptic` | ``E^2 < 1`` | ``R → -∞``: no orbit reaches infinity |
| `:parabolic` | ``E^2 = 1`` | ``R`` is a cubic and grows like ``2r^3`` |
| `:hyperbolic` | ``E^2 > 1`` | ``R → +∞``: orbits come in from infinity |

Outside the outer horizon ``r_+ = 1 + \sqrt{1 - a^2}``, the region ``R ≥ 0`` splits into
intervals bounded by roots of ``R``, by the horizon and by infinity. Each interval is one
possible orbit, and so is each repeated root on which an orbit can sit at constant radius.
[`kerr_geodesic`](@ref) finds all of them and builds one *member* for each. Where the real
roots of ``R`` lie relative to ``r_+``, and how many coincide, decides what the member does;
that arrangement is its *case*.

## The six classes

Every case belongs to one of six classes. A case is named by the class letter and a
number, such as `:A1`, `:K4` or `:C12`.

| Class | Letter | Symbol | Slot of the family | Motion |
| :--- | :---: | :--- | :--- | :--- |
| Stable | A | `:stable` | `Stable` | bound between two turning points, or on a stable circular or spherical orbit |
| Critical | K | `:critical` | `Critical` (a tuple) | on, or asymptotic to, an unstable or marginally stable circular or spherical orbit |
| Plunge | B | `:plunge` | `Plunge` | from a turning point into the black hole |
| Capture | C | `:capture` | `Capture` | from infinity into the black hole, ``E ≥ 1`` |
| Scatter | D | `:scatter` | `Scatter` | from infinity through a turning point and back to infinity, ``E ≥ 1`` |
| Trapped | N | `:trapped` | `Trapped` | ``E < 0``, inside the ergoregion: out of the past horizon, through a turning point, into the future horizon |

The table is available as [`KERR_GEO_CLASSES`](@ref). [`kerr_geo_class`](@ref) returns the
row of a class, a letter or a case, [`kerr_geo_case_class`](@ref) the class of a case, and
[`kerr_geo_member_class`](@ref) the class of a member.

## The cases

In the tables below the real roots of ``R`` are numbered in increasing order, ``r_+`` is the
outer horizon, and a root written twice or three times is a double or triple root.

### Stable (A)

| Case | Energy | Roots of ``R`` | Motion |
| :--- | :--- | :--- | :--- |
| A1 | ``E^2 < 1`` | ``r_1 < r_+ < r_2 < r_3 < r_4`` | oscillates between periapsis ``r_3`` and apoapsis ``r_4`` |
| A2 | ``E^2 < 1`` | ``r_1 < r_+ < r_2 < r_3 = r_4`` | stable circular or spherical orbit on the double root |

### Critical (K)

A member is Critical when its radius sits on, or tends asymptotically to, a repeated root
``r_c > r_+`` that is unstable (a double root with ``R''(r_c) > 0``) or marginally stable (a
triple root). [`kerr_geo_is_critical`](@ref) applies this rule. One such root supports up to
three members, told apart by their *role*, the side of the root on which they lie:

| Roots of ``R`` | Energy | On the root | Outer side | Inner side |
| :--- | :--- | :---: | :---: | :---: |
| ``r_1 < r_+ < r_2 = r_3 = r_4`` | ``E^2 < 1`` | K1 | | K2 |
| ``r_1 < r_+ < r_c = r_c < r_a`` | ``E^2 < 1`` | K3 | K4 | K5 |
| ``r_1 < r_+ < r_2 = r_3`` | ``E^2 = 1`` | K6 | K7 | K8 |
| ``r_1 < r_2 < r_+ < r_3 = r_4`` | ``E^2 > 1`` | K9 | K10 | K11 |

On the root the member is a circular or spherical orbit: the ISCO (equatorial) or ISSO
(inclined) for the triple root K1, and an unstable orbit for K3, K6 and K9, where K6 with
``E = 1`` is the marginally bound orbit. On the outer side, K4 is the homoclinic orbit: it
leaves ``r_c``, turns at the apoapsis ``r_a`` and returns to ``r_c``; K7 and K10 come in from
infinity and settle onto the root. On the inner side, K2, K5, K8 and K11 leave the root inward
and cross the horizon.

The roles are listed in [`CRITICAL_ROLES`](@ref) and returned by
[`kerr_geo_critical_role`](@ref); a member's role is `member.Role`. The role describes a
position, not a direction: K7 and K10 move inward on the outer side. The family keeps its
Critical members in a tuple ordered by role: on the root, outer, inner.

### Plunge (B)

A plunge falls inward from a finite outer turning point. Subextremal members normally put
that point at ``λ = 0`` and reach the future horizon at finite Mino time; axis infall instead
puts the horizon at zero. Exact-extremal crossing members also use a horizon origin, while
B-X1 and B-X2 approach the horizon only asymptotically. See [Where the coordinates are zero](@ref)
for the coordinate origins.

| Case | Energy | Roots of ``R`` | Motion |
| :--- | :--- | :--- | :--- |
| B1 | ``E^2 < 1`` | ``r_1 < r_+ < r_2 < r_3 < r_4`` | from ``r_2`` into the horizon; shares its constants with A1 |
| B2 | ``E^2 < 1`` | ``r_1 < r_+ < r_2 < r_3 = r_4`` | from ``r_2`` into the horizon; shares its constants with A2 |
| B3 | ``E^2 < 1`` | ``r_1 < r_2 < r_3 < r_+ < r_4`` | from ``r_4``, the only root outside the horizon |
| B4 | ``E^2 < 1`` | ``r_1 < r_+ < r_2``, one complex pair | from ``r_2`` into the horizon |
| B5 | ``E^2 = 1`` | ``r_1 < r_+ < r_2 < r_3`` | from ``r_2`` into the horizon; shares its constants with D1 |
| B6 | ``E^2 > 1`` | ``r_1 < r_2 < r_+ < r_3 < r_4`` | from ``r_3`` into the horizon; shares its constants with D2 |
| B7 | ``E^2 < 1`` | ``r_d = r_d < s < r_+ < r_o`` | from ``r_o``; a double root below a simple one inside the horizon |
| B8 | ``E^2 < 1`` | ``s < r_d = r_d < r_+ < r_o`` | from ``r_o``; a simple root below a double one inside the horizon |
| B9 | ``E^2 < 1`` | ``r_d = r_d = r_d < r_+ < r_o`` | from ``r_o``; a triple root inside the horizon, as in radial infall at ``a = 0`` |

### Capture (C)

Every capture comes in from infinity without a turning point and crosses the future
horizon. The cases differ in the roots of ``R``, all of which lie inside the horizon; they
fix the closed form of ``r(λ)``.

| Case | Energy | Roots of ``R`` |
| :--- | :--- | :--- |
| C1 | ``E^2 = 1`` | ``r_1 < r_+``, one complex pair |
| C2 | ``E^2 = 1`` | ``r_1 < r_2 < r_3 < r_+`` |
| C3 | ``E^2 > 1`` | ``r_1 < r_2 < r_+``, one complex pair |
| C4 | ``E^2 > 1`` | ``r_1 < r_2 < r_3 < r_4 < r_+`` |
| C5 | ``E^2 > 1`` | two complex pairs, no real root |
| C6 | ``E^2 = 1`` | ``r_d = r_d < s < r_+`` |
| C7 | ``E^2 = 1`` | ``s < r_d = r_d < r_+`` |
| C8 | ``E^2 = 1`` | ``r_d = r_d = r_d < r_+`` |
| C9 | ``E^2 > 1`` | ``s_1 < r_d = r_d < s_2 < r_+`` |
| C10 | ``E^2 > 1`` | ``s_1 < s_2 < r_d = r_d < r_+`` |
| C11 | ``E^2 > 1`` | ``r_d = r_d < r_+``, one complex pair |
| C12 | ``E^2 > 1`` | ``s < r_d = r_d = r_d < r_+`` |

Orbits from infinity that settle onto an unstable root are Critical (K7, K10).

### Scatter (D)

| Case | Energy | Roots of ``R`` | Motion |
| :--- | :--- | :--- | :--- |
| D1 | ``E^2 = 1`` | ``r_1 < r_+ < r_2 < r_3`` | from infinity, turns at ``r_3``, back to infinity; shares its constants with B5 |
| D2 | ``E^2 > 1`` | ``r_1 < r_2 < r_+ < r_3 < r_4`` | from infinity, turns at ``r_4``, back to infinity; shares its constants with B6 |

### Trapped (N)

With ``E < 0`` a future-directed orbit exists only inside the ergoregion. It leaves the past
(white-hole) horizon, turns at its outer turning point and crosses the future (black-hole)
horizon. The six cases are the six root structures that allow this.

| Case | Energy | Roots of ``R`` | Turning point |
| :--- | :--- | :--- | :--- |
| N1 | ``-1 < E < 0`` | ``r_1 < r_+ < r_2 < r_3 < r_4`` | ``r_2`` |
| N2 | ``-1 < E < 0`` | ``r_1 < r_+ < r_2 < r_3 = r_4`` | ``r_2`` |
| N3 | ``-1 < E < 0`` | ``r_1 < r_2 < r_3 < r_+ < r_4`` | ``r_4``, the only root outside the horizon |
| N4 | ``-1 < E < 0`` | ``r_1 < r_+ < r_2``, one complex pair | ``r_2`` |
| N5 | ``E = -1`` | ``r_1 < r_+ < r_2 < r_3`` (``R`` is a cubic) | ``r_2`` |
| N6 | ``E < -1`` | ``r_1 < r_2 < r_+ < r_3 < r_4`` | ``r_3`` |

[`kerr_energy_sign`](@ref) separates these orbits from the others; the energy regime is set
by ``E^2 - 1`` for either sign of ``E``.

## Constants shared by several members

One set of constants often allows more than one orbit. These go together:

- A1 with B1 and A2 with B2: a stable orbit and the plunge inside it;
- D1 with B5 and D2 with B6: a scattered orbit and the plunge inside it;
- K1 with K2; K3, K4 and K5; K6, K7 and K8; K9, K10 and K11: the Critical members of one
  repeated root.

A family built from such constants contains all of them. `kg.Status.case_ids` lists the
cases, and [`kerr_geo_members`](@ref) returns the members in class order.

## Horizon and extremal tiers

Two degenerate situations produce members outside the numbering above. Their IDs carry the
class letter and a tier letter; in code the IDs use an underscore (`:A_H1`), and
[`kerr_geo_case_name`](@ref) gives the display name (`"A-H1"`).

The **horizon tier** exists for ``0 < |a| < 1`` when ``P(r_+) = E(r_+^2 + a^2) - aL_z = 0``.
Since ``R(r_+) = P(r_+)^2``, the horizon is then itself a root of ``R``.

| Member | Energy | Motion |
| :--- | :--- | :--- |
| A-H1 | ``E < 1`` | oscillates between the two roots outside the horizon |
| A-H2 | ``E < 1`` | stable circular or spherical orbit on a double root outside the horizon |
| D-H1 | ``E = 1`` | from infinity, turns at its only root outside the horizon, back to infinity |
| D-H2 | ``E > 1`` | as D-H1 |

The **extremal tier** exists at ``|a| = 1`` when ``P(r_+) = 2E - aL_z = 0``. The horizon
``r_+ = 1`` is then a double or triple root of ``R``, which an orbit approaches only
asymptotically.

| Member | Energy | Motion |
| :--- | :--- | :--- |
| A-X1 | ``E < 1`` | oscillates between two roots outside the horizon |
| A-X2 | ``E < 1`` | stable spherical orbit |
| B-X1 | ``E < 1`` | turns at its outer root (``λ = 0``) and approaches the double root on the horizon as ``λ → ±∞`` |
| B-X2 | ``E < 1`` | as B-X1, with a triple root on the horizon |
| C-X1 | ``E = 1`` | from infinity toward the double root on the horizon |
| C-X2 | ``E > 1`` | as C-X1 |
| C-X3 | ``E = 1`` | from infinity toward the triple root on the horizon |
| C-X4 | ``E > 1`` | as C-X3 |
| D-X1 | ``E = 1`` | from infinity, turns at its outer root, back to infinity |
| D-X2 | ``E > 1`` | as D-X1 |

The IDs are listed in [`HORIZON_TIER_IDS`](@ref) and [`EXTREMAL_TIER_IDS`](@ref). A member
records its tier in `member.Tier`: `:primary`, `:horizon` or `:extremal`. At ``|a| = 1`` every
member has `Tier == :extremal`, including those with ``P(r_+) ≠ 0``, which keep their
primary case IDs. A repeated root on the horizon is not Critical: the Critical rule asks for
``r_c > r_+``.

## Polar motion

The polar motion is classified separately, from the roots of ``Θ(z)``. The sign of ``Q``
decides; only ``Q = 0`` itself is the equatorial limit, however small a nonzero ``Q`` is. The
constants fix the motion except in two situations, where the keyword `polar_sector` chooses:
with ``Q = 0`` and ``L_z^2 < a^2(E^2 - 1)`` the orbit can lie in the equatorial plane or
approach it asymptotically, and with ``L_z = 0`` and ``Q < 0`` it can pass over the poles or
stay in one hemisphere.

| Sector | Constants | Motion |
| :--- | :--- | :--- |
| `:pendular` | ``Q > 0`` | crosses the equator and turns at the same latitude north and south |
| `:equatorial` | ``Q = 0`` | stays in the equatorial plane |
| `:equator_attractive` | ``Q = 0``, ``L_z^2 < a^2(E^2 - 1)`` | approaches the equatorial plane asymptotically |
| `:vortical` | ``Q < 0``, ``E > 1`` | stays in one hemisphere, oscillating between two latitudes |
| `:constant_latitude` | a double root of ``Θ`` | stays at one latitude off the equator |
| `:axis_crossing` | ``L_z = 0``, ``Q > a^2(1 - E^2)``, ``Q ≠ 0`` | passes over the poles; ``φ`` gains ``π`` at every pass (the limit ``L_z → 0^+``) |
| `:axis_constant` | ``L_z = 0``, ``Q = a^2(1 - E^2)`` | moves along the spin axis |

With ``L_z = Q = 0`` and ``E > 1`` the sector follows the spin. For ``a ≠ 0`` the default is
`:equatorial` and `polar_sector = :axis_crossing` selects motion over the poles (at ``L_z = 0``
the equator-attractive motion is that same axis-crossing motion). For ``a = 0`` the only sector
is `:axis_constant` (radial infall along a fixed direction, case C12), the default with
`axis = :north`. A `polar_sector` the constants do not allow raises an error in both cases.

[`kerr_polar_sector_candidates`](@ref) lists the sectors a set of constants allows, and
`member.Roots.polar.sector` records the one a member uses. Motion confined to one hemisphere
takes `polar_hemisphere = :north` (the default) or `:south`, and motion along the axis takes
`axis = :north` or `:south`.

## The classification itself

[`kerr_geodesic`](@ref) runs the classifier and builds the members. The classifier can also
be called on its own. [`kerr_geo_classify`](@ref) returns a
[`KerrGeoClassification`](@ref) with every admitted region of radial motion, each a
[`KerrGeoRadialComponent`](@ref) with its case, class, endpoints and polar sector;
[`kerr_geo_root_structure`](@ref) returns the roots of ``R`` with their multiplicities; and
[`kerr_geo_case`](@ref) returns the definition of a case as a [`KerrGeoCaseSpec`](@ref).

```@example classes
using KerrGeodesics

a, E, Lz, Q = 0.7, 0.9171300256198305, 2.2591913519439517, 2.898994491984013
classification = kerr_geo_classify(a, E, Lz, Q)
[(c.CaseId, c.BroadClass, c.LowerEndpoint.Radius, c.UpperEndpoint.Radius)
 for c in classification.Components]
```

```@example classes
spec = kerr_geo_case(:K4)
(spec.RootOrdering, spec.AllowedInterval, spec.PastEndpoint, spec.FutureEndpoint)
```

Constants with ``E < 0`` and the horizon and extremal tiers are classified by their own
functions ([`kerr_geo_trapped_classify`](@ref), [`kerr_geo_extremal_family`](@ref)), which
`kerr_geodesic` calls for them.
