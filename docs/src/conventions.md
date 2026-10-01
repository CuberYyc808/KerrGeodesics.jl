# Conventions

## Units and coordinates

Units are ``G = c = M = 1``, with ``M`` the mass of the black hole; the spin satisfies
``|a| ≤ 1``. Coordinates are Boyer–Lindquist ``(t, r, θ, φ)``, with ``z = \cos θ`` used for
the polar motion. The horizons are at ``r_± = 1 ± \sqrt{1 - a^2}``
([`kerr_horizons`](@ref)), and

```math
Σ = r^2 + a^2 z^2, \qquad Δ = r^2 - 2r + a^2 = (r - r_+)(r - r_-).
```

## Mino time and the equations of motion

Every coordinate is a function of Mino time ``λ``, related to proper time by
``dτ/dλ = Σ``. With ``P(r) = E(r^2 + a^2) - aL_z`` ([`kerr_radial_momentum`](@ref)), the
geodesic equations read

```math
\begin{aligned}
\frac{dr}{dλ} &= ±\sqrt{R(r)}, &
\frac{dz}{dλ} &= ±\sqrt{Θ(z)}, \\
\frac{dt}{dλ} &= \frac{(r^2 + a^2)\,P(r)}{Δ} - a^2E\,(1 - z^2) + aL_z, &
\frac{dφ}{dλ} &= \frac{a\,P(r)}{Δ} - aE + \frac{L_z}{1 - z^2},
\end{aligned}
```

with ``R`` and ``Θ`` as in [What decides the motion](@ref). The sign of ``dr/dλ`` is that of
the leg (outgoing or incoming), and ``dθ/dλ = -(dz/dλ)/\sin θ``. Each rate for ``t``, ``φ``
and ``τ`` is a function of ``r`` plus a function of ``z``, so each coordinate is the integral
of a radial part plus the integral of a polar part:

```math
\begin{aligned}
t(λ) &= t(λ_0) + \int_{λ_0}^{λ}\left[\frac{(r^2+a^2)\,P(r)}{Δ} + aL_z - a^2E\,(1 - z^2)\right]dλ', \\
φ(λ) &= φ(λ_0) + \int_{λ_0}^{λ}\left[\frac{a\,P(r)}{Δ} - aE + \frac{L_z}{1 - z^2}\right]dλ', \\
τ(λ) &= τ(λ_0) + \int_{λ_0}^{λ}\bigl[r^2 + a^2 z^2\bigr]\,dλ',
\end{aligned}
```

with ``r = r(λ')``, ``z = z(λ')`` along the orbit and ``λ_0`` the reference event of the member
([Where the coordinates are zero](@ref)). The radial parts are what `radial_t`, `radial_phi`
and `radial_tau` return; over a monotone leg they equal ``\int f(r)\,dr/\sqrt{R(r)}`` with
``f = (r^2+a^2)P/Δ``, ``aP/Δ - aE``, ``r^2`` (see [Radius as the variable](@ref)).

## Constants of motion and APEX parameters

The constants are the energy ``E``, the axial angular momentum ``L_z`` and the Carter
constant ``Q``, all per unit rest mass. They are used exactly as given.

An orbit with a periapsis can also be named by its APEX parameters. The semi-latus rectum
``p`` and the eccentricity ``e`` fix the radial turning points,

```math
r_{\rm apo} = \frac{p}{1 - e}, \qquad r_{\rm peri} = \frac{p}{1 + e},
```

and ``x = \cos ι`` fixes the polar ones: the orbit reaches ``z = ±\sqrt{1 - x^2}``, and the
sign of ``x`` is the sign of ``L_z`` (``x > 0`` for prograde orbits). Bound orbits have
``0 ≤ e < 1``, ``e = 1`` is parabolic and ``e > 1`` hyperbolic. `kerr_geodesic(a, p, e, x)`
converts ``(p, e, x)`` to ``(E, L_z, Q)`` and continues as for constants; a spin within
``8ε`` of ``±1`` is taken as exactly ``±1`` there, and the input is kept in
`Status.input_provenance`.

## Where the coordinates are zero

Each member fixes the zero of ``λ`` at a definite event, and ``τ`` vanishes at ``λ = 0`` for
every member. ``t`` and ``φ`` vanish at ``λ = 0`` too, except where ``λ = 0`` lies on the
future horizon: ``t`` and ``φ`` diverge there, and vanish instead at a reference radius. `ReferenceZero` records all of
this: `lambda0_event`, `t_phi_zero_event` with `t_phi_zero_lambda` and `t_phi_zero_radius`,
`tau_zero_event`, and `lambda_regular`, the ``λ`` where ``v`` and ``ψ`` vanish.

| Members | ``λ = 0`` | ``t = φ = 0`` | ``v = ψ = 0`` |
| :--- | :--- | :--- | :--- |
| A1, A2 | periapsis and the northern polar turning point | ``λ = 0`` | |
| A-H1, A-X1 | the inner turning point | ``λ = 0`` | |
| K1, K3, K6, K9, A-H2, A-X2 | the polar reference event | ``λ = 0`` | |
| K4 | apoapsis | ``λ = 0`` | |
| K7, K10 | the reference radius | ``λ = 0`` | ``λ = 0`` |
| K2, K5, K8, K11 | the future horizon | the reference radius | ``λ = 0`` |
| B1–B8 | the turning point | ``λ = 0`` | the future horizon |
| B9 and motion along the axis | the future horizon | the reference radius | ``λ = 0`` |
| C1–C12 | the future horizon | the reference radius | ``λ = 0`` |
| D1, D2, D-H1, D-H2 | the turning point | ``λ = 0`` | |
| N1–N6 | the turning point | ``λ = 0`` | the future horizon (``u = χ = 0`` on the past horizon) |
| B-X1, B-X2, D-X1, D-X2 | the turning point | ``λ = 0`` | |
| C-X1…C-X4 | the reference radius | ``λ = 0`` | |

`reference_radius` sets the reference radius. Its default is ``λ_∞/2``, halfway in Mino time
between infinity and the horizon, for Capture members; ``r_+ + 0.55\,(r_c - r_+)`` for K2,
K8 and K11 and ``(r_+ + r_c)/2`` for K5; and ``r_c + \max(1, r_c - r_+)`` for K7 and K10.
The initial phases of a Stable member, `initPhases = (qt0, qr0, qθ0, qφ0)`, shift its
``t``, radial phase, polar phase and ``φ`` at ``λ = 0``.

## The polar phase

The keyword `polar_phase` sets the phase of the polar motion at the member's reference
event, and `member.ReferenceZero.polar_phase_convention` says what phase 0 means:

| Convention | Phase 0 |
| :--- | :--- |
| `:northern_turning_point_at_zero_phase` | the northern turning point, ``z = +z_{max}`` |
| `:equator_northward_at_zero_phase` | the equator, crossed northward |
| `:maximum_absolute_latitude_at_zero_phase` | the turning point farthest from the equator (vortical motion) |
| `:not_applicable` | no polar oscillation: equatorial motion or motion along the axis |

For K2, K8 and K11 the phase is given at the reference-radius event
(`ReferenceZero.polar_phase_event = :reference_radius`) rather than at ``λ = 0``, which lies
on the horizon.

## Coordinates regular at the horizon

``t`` and ``φ`` diverge logarithmically on a horizon. With the tortoise coordinate

```math
r_* = r + \frac{2r_+}{r_+ - r_-}\ln\frac{r - r_+}{2} - \frac{2r_-}{r_+ - r_-}\ln\frac{r - r_-}{2},
\qquad \frac{dr_*}{dr} = \frac{r^2 + a^2}{Δ},
```

([`kerr_rstar`](@ref)) and the horizon azimuth

```math
φ_H = \frac{a}{r_+ - r_-}\ln\left|\frac{r - r_+}{r - r_-}\right|, \qquad \frac{dφ_H}{dr} = \frac{a}{Δ},
```

the ingoing coordinates ``v = t + r_*`` and ``ψ = φ + φ_H`` are finite on the future
horizon, and the outgoing coordinates ``u = t - r_*`` and ``χ = φ - φ_H`` are finite on the
past horizon. At ``|a| = 1`` both horizons sit at ``r = 1`` and ``r_*`` takes its limiting
form ``r + 2\ln(r - 1) - 2/(r - 1) - 2\ln 2``.

A member that crosses the future horizon returns ``v``, ``ψ`` there directly; they are
computed without the cancellation between ``t`` and ``r_*``, so they keep full precision next
to the horizon.

## The sign of the spin

The reflection ``(a, L_z, φ) → (-a, -L_z, -φ)``, at fixed ``E``, ``Q``, Mino origin and polar
phase, maps a member at spin ``a`` onto the member at ``-a``: ``r``, ``z``, ``t``, ``τ``,
``r_*``, ``u`` and ``v`` are unchanged and ``φ``, ``ψ`` and ``χ`` change sign.

## Energy

[`kerr_energy_regime`](@ref) returns the sign of ``E^2 - 1`` as `:elliptic`, `:parabolic` or
`:hyperbolic`, for either sign of ``E``; [`kerr_energy_sign`](@ref) returns the sign of ``E``
itself. Future-directed motion with ``E < 0`` exists only inside the ergoregion, which is the
Trapped class.
