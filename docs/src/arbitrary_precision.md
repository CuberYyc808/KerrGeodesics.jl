# Arbitrary precision

By default every computation is in `Float64`. The package also computes in `BigFloat`, at any
number of bits: the precision is set by the input, and nothing else changes in how an orbit is
built or used.

## Usage

There are three ways to ask for more digits.

**`BigFloat` input.** Every number is computed in the floating-point type of the input.
`BigFloat` constants give a `BigFloat` orbit at their precision:

```julia
setprecision(BigFloat, 256) do
    kerr_geodesic(big"0.9", (big"0.94", big(2), big(4)))
end
```

**The `precision` keyword.** `precision = p` converts the input to `p`-bit `BigFloat`:

```@example precision
using KerrGeodesics

kg64 = kerr_geodesic(0.9, (0.94, 0.1, 12.0))                    # Float64, the default
kg = kerr_geodesic(0.9, (0.94, 0.1, 12.0); precision=256)
λ = 0.5
kg.Plunge.Trajectory.r(big(λ))
```

```@example precision
kg.Plunge.Trajectory.r(big(λ)) - kg64.Plunge.Trajectory.r(λ)     # the Float64 error
```

**Rational input.** With `precision`, a rational number is rounded once, in the target
precision, so `9//10` is ``9/10`` to the last of the 256 bits:

```@example precision
kgq = kerr_geodesic(9//10, (47//50, 1//10, 12); precision=256)
kgq.Plunge.ConstantsOfMotion.E
```

The same keyword is accepted by `kerr_geodesic(a, p, e, x)`, `kerr_geo_orbit`,
`kerr_geo_stable`, `kerr_geo_plunge`, `kerr_geo_frequencies` and
`kerr_geo_constants_of_motion`. Every other function follows the type of its arguments.

```@example precision
f = kerr_geo_frequencies(0.9, 10.0, 0.5, 0.8; precision=128)
f["ϒr"]
```

## An orbit remembers its precision

The functions of an orbit evaluate at the precision the orbit was built with, whatever
`setprecision` is in effect where they are called:

```@example precision
r = setprecision(BigFloat, 64) do
    kg.Plunge.Trajectory.r(big(λ))
end
precision(r)
```

Close to the separatrix `Float64` loses digits in the radial frequencies; `BigFloat` restores
them in proportion to the number of bits (see [Numerics and accuracy](accuracy.md)).

## What follows the precision

Every step is computed in the type of the input: the roots of the radial potential and the
classification, the elliptic integrals and Jacobi functions (implemented in the package), the
series, the Chebyshev tables of ``t``, ``φ`` and ``τ`` (their order grows with the number of
digits) and the tolerances, which keep their place between the rounding and the physical scale.
At 256 bits the elliptic functions and the radial roots are accurate to about ``10^{-76}``, the
orbits to about ``10^{-70}``.

## Three points about the input

- `BigFloat(0.9)` is the binary value of the `Float64` number `0.9`, not ``9/10``. To give a
  decimal value, pass a rational (`9//10`) with `precision`, or the string `big"0.9"`.
- `big"0.9"` is parsed when the code is compiled, at the precision in effect then. Inside a
  function, use `parse(BigFloat, "0.9")` or a rational.
- A repeated root or a root on the horizon is a property of the exact constants, recognised to
  within one unit in the last place of the input. Constants rounded to `Float64` from an exact
  degenerate set are degenerate in `Float64` but generally not at 256 bits: compute such
  constants in the target precision.

## Cost

A `BigFloat` orbit is about a thousand times slower than a `Float64` one. For a plunge at
256 bits, building the orbit takes one to two seconds and an evaluation of ``t`` or ``r`` about
0.1–1 ms (`Float64`: about a millisecond and a microsecond); the cost grows roughly with the
number of bits.

## Supported types

`Float64` (the default) and `BigFloat`. Other formats (double-double, quadruple, multi-word
floats) are not supported: the exact zero tests and error-free products of the classification
rely on a binary, correctly rounded type with a full exponent range and an exact `fma`.
