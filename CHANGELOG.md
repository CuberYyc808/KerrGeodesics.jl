# Changelog

## 0.5.0

### New
- **Arbitrary precision.** Every computation follows the floating-point type of the input;
  `BigFloat` input, or the `precision` keyword of `kerr_geodesic`, `kerr_geo_orbit`,
  `kerr_geo_stable`, `kerr_geo_plunge`, `kerr_geo_frequencies` and
  `kerr_geo_constants_of_motion`, gives orbits at any number of bits. Rational input is rounded
  once, in the target precision. An orbit evaluates at the precision it was built with. See
  *Arbitrary precision* in the documentation.
- **APEX turning points.** For a bound eccentric orbit given by `(p, e, x)` the Stable member
  keeps the turning points `p/(1 ∓ e)` and takes the inner roots from the quadratic factor of the
  radial potential of the same constants. Near the separatrix the rounded constants alone can
  merge `r2` and `r3`, which left such orbits without a Stable member (classified K3/K4/K5);
  these orbits now have their Stable member, and the error of the radial frequency grows as
  `ε/(p − p_s)` instead of `ε/(p − p_s)²`. `Status.apex_root_geometry` records whether the
  turning points were used and, if not, why; `Status.component_root_models` gives the root model
  of every component.

### Changed
- **Orbits over the spin axis (`x = 0`, `Lz = 0`).** The azimuth is the limit `Lz → 0⁺`
  (`x → 0⁺`): `φ` gains `π` at every pass over the axis, and the azimuthal frequency `ϒφ` (and
  `Ωφ`) includes `ϒθ`. Before, `φ` had no jump at the axis and `ϒφ` lacked `ϒθ`, which placed
  the orbit after a pass on the mirrored meridian. Applies to the members, `kerr_geo_frequencies`
  and `kerr_geo_orbit`; `−0.0` gives the limit `x → 0⁻`.
- Classification and horizon coordinates: repeated roots are read from the constants within one
  unit of input rounding, and interior clusters by their effect on `R(r₊)`.

### Fixed
- `BigFloat` orbits with `0 < |Lz| ≲ 10⁻¹⁸` (near the axis) failed to build their Stable and
  Capture members at 128 bits and above; the Chebyshev bisection depth now grows with the
  precision (60 in Float64, as before).

### Migration
- `ϒφ`, `Ωφ`, `φ(λ)` of orbits with `x = 0` change as described above; for `x → 0⁻` pass `−0.0`.
- APEX orbits near the separatrix may now have a Stable member together with the critical
  members of the rounded constants; select with `case_id` as before.
- Float64 results of Stable members built from APEX input change at the level of rounding
  (their inner roots now come from the turning-point factorization); constants input is
  unchanged.

### Known issues
- At `|a| = 1`, an interior repeated radial root (inside the horizon) whose constants are a few
  units of rounding away from an exact repeated root can be read as two simple roots: for
  `(a, E) = (±1, 1.2)` with a double root at `r = 0.1`, constants computed in Float64 or at
  128 bits give C4 instead of C10; computed at 256 bits they give C10.
- Near the separatrix the Float64 accuracy of the radial frequency of APEX orbits is limited by
  the rounding of the constants (see *Numerics and accuracy*): its relative error reaches about
  `2 × 10⁻⁷` at `p − p_s = 10⁻⁹` and `10⁻⁶` at `p − p_s = 10⁻¹⁰`. Use `BigFloat` closer to the
  separatrix.
