# Examples

## Stable Bound Orbit

Stable bound geodesics use APEX-like parameters `(a,p,e,x)`.

```julia
using KerrGeodesics

geo = kerr_geodesic(0.9, 10.0, 0.5, 0.8)
stable = geo.Stable

t = stable.Trajectory.t
r = stable.Trajectory.r
theta = stable.Trajectory.θ
phi = stable.Trajectory.ϕ
```

The direct stable constructor is:

```julia
stable = kerr_geo_stable(0.9, 10.0, 0.5, 0.8;
                         initPhases=(0.0, 0.0, 0.0, 0.0))
```

Stable phase offsets are `(qt0, qr0, qtheta0, qphi0)`.

## Finite-Start Plunge From Constants

Finite-start plunge trajectories use constants of motion `(a,E,Lz,Q)`.

```julia
plunge = kerr_geo_plunge(0.9, 0.94, 0.1, 12.0;
                         radial_start=:turning_point)

r0 = plunge.Trajectory.r(0.0)
u0 = plunge.Trajectory.u(0.0)
v0 = plunge.Trajectory.v(0.0)
rstar0 = plunge.Trajectory.rstar(0.0)
```

The unified constructor returns a family object:

```julia
family = kerr_geodesic(0.9, (0.94, 0.1, 12.0);
                       radial_start=:turning_point)
family.RootClass
family.Plunge.Status
```

## Initial Position Controls

Plunge trajectories can be initialized with explicit Mino-time phases:

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

For equatorial or otherwise degenerate polar sectors, `initial_theta` defaults
to `pi/2` and the polar phase is set to zero.

## Horizon-Aligned Advanced Time

For finite-start plunge waveform diagnostics, the time origin can be shifted so
the future-horizon advanced-time anchor satisfies $v_H=0$:

```julia
plunge = kerr_geo_plunge(a, E, Lz, Q;
                         time_origin=:future_horizon_v_zero)

plunge.Status.horizon_time_shift
plunge.Status.horizon_v_anchor_method
```

This shifts the additive coordinate-time origin. It does not make the retarded
time finite at the future horizon; $u=v-2r_*$ diverges as $r_*\to-\infty$.

When available, the near-horizon series callables evaluate the same convention
as functions of `rstar`:

```julia
lambda_probe = 10.0
rstar_cutoff = plunge.Trajectory.rstar(lambda_probe)
v_series = plunge.Trajectory.v_rstar_series(rstar_cutoff)
u_series = plunge.Trajectory.u_rstar_series(rstar_cutoff)
```

The series evaluator is a near-horizon bridge at fixed `rstar`; it does not
replace the trajectory path from `lambda` to `r` and `rstar`.

## Root Classification

The plunge branch label is available as:

```julia
plunge.OrbitClass
family.RootClass
```

Branch labels are implementation metadata and should be propagated into
downstream trajectory, source, and waveform manifests.

For stable-orbit classification, use:

```julia
metadata = kerr_geo_orbit_type_metadata(a, p, e, x)
metadata.labels
metadata.at_separatrix
```

Inputs inside the current separatrix roundoff guard are evaluated at the
separatrix radius and labeled `Separatrix`.
