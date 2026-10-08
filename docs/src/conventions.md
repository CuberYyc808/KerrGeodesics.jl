# Conventions

## Units and coordinates

Units are ``G = c = M = 1``, with ``M`` the mass of the black hole; the spin satisfies
``|a| \leq 1``. Coordinates are Boyer–Lindquist ``(t, r, \theta, \phi)``, with ``z = \cos \theta`` used for
the polar motion. The horizons are at ``r_\pm = 1 \pm \sqrt{1 - a^2}``
([`kerr_horizons`](@ref)), and

```math
\Sigma = r^2 + a^2 z^2, \qquad \Delta = r^2 - 2r + a^2 = (r - r_+)(r - r_-).
```

## Mino time and the equations of motion

Every coordinate is a function of Mino time ``\lambda``, related to proper time by
``d\tau/d\lambda = \Sigma``. With ``P(r) = E(r^2 + a^2) - aL_z`` ([`kerr_radial_momentum`](@ref)), the
geodesic equations read

```math
\begin{aligned}
\frac{dr}{d\lambda} &= \pm\sqrt{R(r)}, &
\frac{dz}{d\lambda} &= \pm\sqrt{\Theta(z)}, \\
\frac{dt}{d\lambda} &= \frac{(r^2 + a^2)\,P(r)}{\Delta} - a^2E\,(1 - z^2) + aL_z, &
\frac{d\phi}{d\lambda} &= \frac{a\,P(r)}{\Delta} - aE + \frac{L_z}{1 - z^2},
\end{aligned}
```

with ``R`` and ``\Theta`` as in [What decides the motion](@ref). The sign of ``dr/d\lambda`` is that of
the leg (outgoing or incoming), and ``d\theta/d\lambda = -(dz/d\lambda)/\sin \theta``. Each rate for ``t``, ``\phi``
and ``\tau`` is a function of ``r`` plus a function of ``z``, so each coordinate is the integral
of a radial part plus the integral of a polar part:

```math
\begin{aligned}
t(\lambda) &= t(\lambda_0) + \int_{\lambda_0}^{\lambda}\left[\frac{(r^2+a^2)\,P(r)}{\Delta} + aL_z - a^2E\,(1 - z^2)\right]d\lambda', \\
\phi(\lambda) &= \phi(\lambda_0) + \int_{\lambda_0}^{\lambda}\left[\frac{a\,P(r)}{\Delta} - aE + \frac{L_z}{1 - z^2}\right]d\lambda', \\
\tau(\lambda) &= \tau(\lambda_0) + \int_{\lambda_0}^{\lambda}\bigl[r^2 + a^2 z^2\bigr]\,d\lambda',
\end{aligned}
```

with ``r = r(\lambda')``, ``z = z(\lambda')`` along the orbit and ``\lambda_0`` the reference event of the member
([Where the coordinates are zero](@ref)). The radial parts are what `radial_t`, `radial_phi`
and `radial_tau` return; over a monotone leg they equal ``\int f(r)\,dr/\sqrt{R(r)}`` with
``f = (r^2+a^2)P/\Delta``, ``aP/\Delta - aE``, ``r^2`` (see [Radius as the variable](@ref)).

At ``L_z = 0`` the orbit passes over the spin axis, where ``\phi`` is undefined. ``\phi`` is the
limit ``L_z \to 0^+`` (``x \to 0^+``): it gains ``\pi`` at every pass over the axis, so the Cartesian
position continues to the opposite meridian, and the azimuthal frequency ``\Upsilon_\phi`` includes
``\Upsilon_\theta``.

## Constants of motion and APEX parameters

The constants are the energy ``E``, the axial angular momentum ``L_z`` and the Carter
constant ``Q``, all per unit rest mass. They are used exactly as given.

An orbit with a periapsis can also be named by its APEX parameters. For a bound orbit,
the semi-latus rectum ``p`` and eccentricity ``0 \leq e < 1`` fix the radial turning points,

```math
r_{\rm apo} = \frac{p}{1 - e}, \qquad r_{\rm peri} = \frac{p}{1 + e},
```

whereas ``r_{\rm peri} = p/(1+e)`` also applies to scattering orbits. At ``e = 1`` the
outer endpoint is infinity; for ``e > 1``, ``p/(1-e)`` is not a physical apoapsis.
The inclination parameter ``x = \cos \iota`` fixes the polar turning points
``z = \pm\sqrt{1 - x^2}``, and its sign is the sign of ``L_z``. For nonzero spin, motion is
prograde when ``aL_z > 0`` (equivalently ``ax > 0``), and retrograde when ``aL_z < 0``.
At ``a = 0`` these labels have no distinction relative to the black-hole spin.
``e = 1`` is parabolic and ``e > 1`` hyperbolic. `kerr_geodesic(a, p, e, x)`
converts ``(p, e, x)`` to ``(E, L_z, Q)`` and continues as for constants; a spin within
``8\epsilon`` of ``\pm1`` is taken as exactly ``\pm1`` there, and the input is kept in
`Status.input_provenance`.

For a bound eccentric orbit (``|a| < 1``, ``0 < e < 1``, ``0 < E < 1``) the turning points
``p/(1 \mp e)`` are part of the input, and the Stable member keeps them: its roots are
``r_1 = p/(1-e)``, ``r_2 = p/(1+e)`` and the two roots of the quadratic factor
``R(r)/\bigl((r - r_1)(r - r_2)\bigr)`` of the same constants. Near the separatrix the
constants alone, rounded to the working precision, can merge ``r_2`` and ``r_3`` into a
repeated root, and the bound orbit would be classified as critical. The other members use
the roots of the constants. `Status.apex_root_geometry` records whether the turning points
were used (`accepted`) or, if not, why (`reason`, e.g. `:gap_unresolved` when ``r_2 - r_3``
is within rounding), and `Status.component_root_models` gives the root model of each
component.

## Where the coordinates are zero

Each member fixes the zero of ``\lambda`` at a definite event, and ``\tau`` vanishes at ``\lambda = 0``.
With the default phases, ``t`` and ``\phi`` usually vanish there too. Subextremal members
anchored at the future horizon instead fix their Boyer–Lindquist zeros at an exterior
reference radius. Exact-extremal crossing members fix their constants through the regular
ingoing chart, rather than through a finite zero of ``t`` and ``\phi``.
`ReferenceZero` records the event and coordinate conventions in `lambda0_event`,
`t_phi_zero_event`, `t_phi_zero_lambda`, `t_phi_zero_radius`, `tau_zero_event` and
`lambda_regular`. The last field gives the zero of ``v`` and ``\psi`` when it is defined.
The table assumes default phases and `phi0 = 0`; motion along the axis with nonzero
`phi0` keeps that constant azimuth.

| Members | ``\lambda = 0`` | ``t = \phi = 0`` | ``v = \psi = 0`` |
| :--- | :--- | :--- | :--- |
| A1, A2 | periapsis and the northern polar turning point | ``\lambda = 0`` | |
| A-H1, A-X1 | the inner turning point | ``\lambda = 0`` | |
| K1, K3, K6, K9, A-H2, A-X2 | the polar reference event | ``\lambda = 0`` | |
| K4 | apoapsis | ``\lambda = 0`` | |
| K7, K10 | the reference radius | ``\lambda = 0`` | ``\lambda = 0`` |
| K2, K5, K8, K11, ``\lvert a\rvert < 1`` | the future horizon | the reference radius | ``\lambda = 0`` |
| B1 to B9, ``\lvert a\rvert < 1``, except axis infall | the turning point | ``\lambda = 0`` | the future horizon |
| Axis infall, including Schwarzschild B9 | the future horizon | the reference radius for ``t``; constant ``\phi`` | ``\lambda = 0`` |
| C1 to C12, ``\lvert a\rvert < 1``, except axis infall | the future horizon | the reference radius | ``\lambda = 0`` |
| Primary Plunge, Capture and inward Critical crossing members, ``\lvert a\rvert = 1``, ``P(r_+) \neq 0`` | the future horizon | fixed through the ingoing chart | ``\lambda = 0`` |
| D1, D2, D-H1, D-H2 | the turning point | ``\lambda = 0`` | |
| N1–N6 | the turning point | ``\lambda = 0`` | the future horizon (``u = \chi = 0`` on the past horizon) |
| B-X1, B-X2, D-X1, D-X2 | the turning point | ``\lambda = 0`` | |
| C-X1…C-X4 | the reference radius | ``\lambda = 0`` | |

For the exact-extremal crossing members in the table,
`t_phi_zero_event = :regular_chart_at_future_horizon`: ``t = v - r_*`` and
``\phi = \psi - \phi_H`` inherit their constants from ``v = \psi = 0`` on the horizon. Their
`t_phi_zero_lambda` and `t_phi_zero_radius` are NaN because no finite zero is prescribed.

`reference_radius` sets the exterior reference radius where that convention is used.
For subextremal Capture members its default is ``r(\lambda_\infty/2)``, where ``\lambda_\infty`` is the Mino-time
endpoint at infinity and ``\lambda = 0`` is the horizon: the reference event lies halfway between
them in Mino time. The defaults are ``r_+ + 0.55\,(r_c - r_+)`` for K2,
K8 and K11 and ``(r_+ + r_c)/2`` for K5; and ``r_c + \max(1, r_c - r_+)`` for K7 and K10.
The initial phases of a Stable member, `initPhases = (qt0, qr0, qθ0, qφ0)`, shift its
``t``, radial phase, polar phase and ``\phi`` at ``\lambda = 0``.

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
(`ReferenceZero.polar_phase_event = :reference_radius`) rather than at ``\lambda = 0``, which lies
on the horizon.

## Coordinates regular at the horizon

At a subextremal horizon, ``t`` has a logarithmic divergence on a crossing trajectory;
the azimuth generally does too, although its divergent coefficient can vanish. For
``|a| < 1`` the tortoise coordinate is

```math
r_* = r + \frac{2r_+}{r_+ - r_-}\ln\frac{r - r_+}{2} - \frac{2r_-}{r_+ - r_-}\ln\frac{r - r_-}{2},
\qquad \frac{dr_*}{dr} = \frac{r^2 + a^2}{\Delta},
```

([`kerr_rstar`](@ref)), and the horizon azimuth is

```math
\phi_H = \frac{a}{r_+ - r_-}\ln\left|\frac{r - r_+}{r - r_-}\right|, \qquad \frac{d\phi_H}{dr} = \frac{a}{\Delta},
```

The ingoing coordinates ``v = t + r_*`` and ``\psi = \phi + \phi_H`` are finite on the future
horizon, and the outgoing coordinates ``u = t - r_*`` and ``\chi = \phi - \phi_H`` are finite on the
past horizon. At ``|a| = 1`` both horizons sit at ``r = 1``. Outside the horizon the
limiting forms, with the package's additive constants, are

```math
r_* = r + 2\ln(r - 1) - \frac{2}{r - 1} - 2\ln 2,
\qquad \phi_H = -\frac{a}{r - 1}.
```

Thus generic exact-extremal crossing has a pole as well as a logarithm in ``t``, and a
pole in ``\phi``; the subextremal logarithmic statement does not apply there.

A member that crosses the future horizon returns ``v``, ``\psi`` there directly; they are
computed without subtracting the divergent Boyer–Lindquist terms. Their accuracy does not
follow from the accuracy of ``t`` and ``\phi`` separately (see [Limits of double precision](@ref)).

## The sign of the spin

The reflection ``(a, L_z, \phi) \to (-a, -L_z, -\phi)``, at fixed ``E``, ``Q``, Mino origin and polar
phase, maps a member at spin ``a`` onto the member at ``-a``: ``r``, ``z``, ``t``, ``\tau``,
``r_*``, ``u`` and ``v`` are unchanged and ``\phi``, ``\psi`` and ``\chi`` change sign.

## Energy

[`kerr_energy_regime`](@ref) returns the sign of ``E^2 - 1`` as `:elliptic`, `:parabolic` or
`:hyperbolic`, for either sign of ``E``; [`kerr_energy_sign`](@ref) returns the sign of ``E``
itself. Future-directed motion with ``E < 0`` exists only inside the ergoregion, which is the
Trapped class.
