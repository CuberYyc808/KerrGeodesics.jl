# Numerics and accuracy

## How the coordinates are computed

``r(λ)`` and ``z(λ)`` are closed forms: Jacobi elliptic functions or elementary functions,
one formula per root structure (`member.Status.formula_family` names it). The formulas
reduce cancellation near turning points and endpoints. Their accuracy still depends on
the constants, phase and intermediate calculations; there is no uniform relative-error
bound over all parameter limits.

``t``, ``φ``, ``τ``, ``v`` and ``ψ`` are integrals over ``λ`` of the rates in
[Mino time and the equations of motion](@ref). Each rate is a radial part plus a polar
part, and each part is tabulated once as a piecewise Chebyshev series fitted to machine
precision. The pieces that diverge or grow without bound are added in closed form, so the
series only represent bounded remainders: the singular parts of ``r_*`` and ``φ_H`` at a
horizon, the growth of ``t`` and ``φ`` towards infinity, and the ``L_z/(1 - z^2)`` spike of
``dφ/dλ`` on orbits that pass close to the spin axis. For periodic polar motion the series
cover one quarter period and the rest follows by symmetry; nonperiodic polar limits are
handled separately.

The polar tables, and the radial tables of Stable members, are built with the member. The
other radial tables are built the first time ``t``, ``φ``, ``τ``, ``v`` or ``ψ`` is
evaluated. At ``|a| = 1`` the radial parts come from closed-form Mino-time integrals of the
radial models instead.

## What a member reports

`member.Status.spectral` describes the member's Chebyshev tables: `achieved` is the largest
error estimate relative to the local size of the tabulated rate, and `pieces` the number of
pieces. A table that cannot reach its tolerance raises an error when it is built, so
`unresolved` is always 0.

```@example accuracy
using KerrGeodesics

stable = kerr_geodesic(0.9, 10.0, 0.5, 0.8).Stable
(achieved = stable.Status.spectral.achieved, pieces = stable.Status.spectral.pieces)
```

A Stable member also reports, in `Status.precision`, how much ``t`` advances over one radial
period and the spacing of double-precision numbers at that value:

```@example accuracy
stable.Status.precision
```

## Arbitrary precision

Everything above holds in the floating-point type of the input. In `BigFloat` the tolerances,
the series and the Chebyshev orders follow the number of digits; see
[Arbitrary precision](arbitrary_precision.md).

## Limits of double precision

Accuracy can be limited both by the conditioning of Float64 inputs and by rounding in
intermediate calculations. The spectral estimate measures the fitted remainder, not the
total coordinate error. The following situations distinguish these effects.

**Long times in the polar phase.** ``z(λ)`` depends on the phase ``u_0 + ωλ``, and rounding
this phase leaves an absolute error of about ``ε|ωλ|``, with ``ε`` the machine epsilon. The
error reaches order one at ``|λ| ∼ 1/(ωε)``.

**Nearly parabolic stable orbits.** As ``E → 1`` the radial period grows without bound. After
one period ``t`` is known only to `Status.precision.t_ulp_per_period`, about ``10^3`` for
``|E - 1| ∼ 10^{-13}``.

**Near the spin axis.** An orbit with small ``L_z`` swings through ``φ`` in a narrow spike
each time it passes the axis. On the spike, ``φ(λ)`` carries the error
``|dφ/dλ|\,\mathrm{ulp}(λ)`` from the rounding of ``λ`` (about ``10^{-3}`` for
``L_z = 10^{-12}``). The change of ``φ`` over a polar period, measured from any point off
the spike, agrees to about ``10^{-15}`` in the tested examples; this is not a bound for
every choice of constants and phase.

**Far from the hole.** For ``E > 1``, ``r(λ) ≈ 1/\bigl(\sqrt{E^2 - 1}\,|λ - λ_∞|\bigr)``
near the end ``λ_∞`` of the Mino domain, so ``λ`` resolves radii only up to about
``1/\bigl(\sqrt{E^2 - 1}\,\mathrm{ulp}(λ_∞)\bigr) ∼ 10^{16}``. The radius-based functions
`radial_*_increment(r1, r2)` take radii directly and avoid this particular resolution limit.

**Boyer–Lindquist time very close to the horizon.** Rounding the intermediate radius can
be amplified in a divergent coordinate even when ``λ`` itself is well resolved. For the
K8 member with ``(a,E,L_z,Q) = (0.7,1,-0.7,16)``, comparison with independent high-precision
quadrature gives relative discrepancies in ``t`` of about ``8 × 10^{-10}`` at
``λ = -10^{-8}`` and ``5 × 10^{-8}`` at ``λ = -10^{-10}``, with the horizon at zero.
This concerns the divergent Boyer–Lindquist coordinate, not an error bound for the regular
coordinates ``v`` and ``ψ``.

**Near the separatrix.** As ``p`` approaches the separatrix ``p_s``, the radial roots ``r_2``
and ``r_3`` merge and the frequencies become sensitive to the constants. For APEX input the
Stable member keeps the turning points (see [Conventions](conventions.md)), so the relative
error of ``ϒ_r`` grows about as ``ε/(p - p_s)``; for a tuple of constants it grows about as
``ε/(p - p_s)^2``, because the rounding of the constants also moves ``r_2``. In Float64,
comparison with independent 120-digit references for 436 orbits (``|a| ≤ 0.99``,
``e = 0.1, 0.5, 0.8``, all inclinations) gives a largest relative error of ``ϒ_r`` of
``1.4 × 10^{-10}``, ``1.2 × 10^{-9}``, ``1.9 × 10^{-8}``, ``2.2 × 10^{-7}`` and
``1.1 × 10^{-6}`` at ``p - p_s = 10^{-6}, 10^{-7}, 10^{-8}, 10^{-9}, 10^{-10}``; the largest
errors are at ``a = 0.99``, where the constants themselves carry a few units of rounding.
Closer to the separatrix, or for tighter requirements, use `BigFloat`.

**Very small nonzero spin.** Some combinations approaching the Schwarzschild limit lose
precision through cancellation. This is an implementation limitation in addition to
input conditioning, and is not described by the spectral fit estimate alone.

## Checking an orbit

[`kerr_geo_diagnose`](@ref) checks a member, or every member of a family, against the
geodesic equations using only its own functions. On samples across the domain it compares
finite differences of the trajectory with ``R``, ``Θ`` and the rates of ``t`` and ``φ``, and
over one polar period it compares the change of ``φ`` with an independent quadrature. It
also times single evaluations.

```@example accuracy
report = kerr_geo_diagnose(kerr_geodesic(0.9, 10.0, 0.5, 0.8))
all(r.ok for r in report)
```

```@example accuracy
[(r.member, round.((r.radial, r.polar, r.time, r.azimuth); sigdigits = 2)) for r in report]
```

The discrepancies are relative and include finite-difference error as well as any error
in the trajectory. Each record also gives the allowance granted for finite differences, the
samples that a finite difference cannot decide (`skipped`), and the timings.

## Speed

After compilation for the relevant member types, building a family, classification
included, takes about 0.1 ms in the measured examples. The first evaluation of
``t`` or ``φ`` builds the radial tables, which takes a few hundredths of a millisecond, up to
about 0.4 ms for orbits that reach infinity. After that a coordinate costs a fraction of a
microsecond per evaluation. These are warmed timings: loading and the first specialization
of a new member type add compilation cost. [`kerr_geo_sample`](@ref) evaluates the nine quantities of the
trajectory and velocity at one ``λ`` in about 2 µs, without the dispatch cost of calling
each function from untyped code.
