# API reference

Every exported name, grouped by purpose. [Working with an orbit](@ref) describes the fields
of the family and member objects in context.

## Building orbits

```@docs
KerrGeodesics
kerr_geodesic
KerrGeodesicFamily
kerr_geo_members
KerrGeoComponent
kerr_geo_member_class
kerr_geo_sample
```

`KerrGeoStableComponent`, `KerrGeoCriticalComponent`, `KerrGeoPlungeComponent`,
`KerrGeoCaptureComponent`, `KerrGeoScatterComponent` and `KerrGeoTrappedComponent` are the
six classes of [`KerrGeoComponent`](@ref), `KerrGeoComponent{:stable}` and so on.

## Members by class

The class constructors `kerr_geo_*_component` and [`kerr_geo_trapped`](@ref) take
`(a, E, Lz, Q)`, classify the constants and return one member. They also take
`(a, (E, Lz, Q))`. When the constants allow more than one member of the class, `case_id`
chooses. Specialized constructors use the independent constants of their sector:
axis infall takes `(a, E; axis=...)`, and equator-attractive capture takes `(a, E, Lz)`
with ``Q = 0``. Their docstrings give the individual signatures and keyword defaults.

### Stable

```@docs
kerr_geo_stable_component
kerr_geo_stability_metadata
```

### Critical

```@docs
kerr_geo_critical_component
kerr_geo_critical_spherical
kerr_geo_critical_homoclinic
kerr_geo_critical_plunge
```

### Plunge

```@docs
kerr_geo_plunge_component
kerr_geo_plunge_axis_infall
```

### Capture

```@docs
kerr_geo_capture_component
kerr_geo_capture_four_complex
kerr_geo_capture_vortical
kerr_geo_capture_constant_latitude
kerr_geo_capture_equator_attractive
kerr_geo_capture_axis_infall
```

### Scatter

```@docs
kerr_geo_scatter_component
kerr_geo_scatter_asymptotic_diagnostics
kerr_geo_scatter_asymptotic_state
```

### Trapped

```@docs
kerr_geo_trapped
kerr_geo_trapped_case
kerr_geo_trapped_classify
KerrGeoTrappedClassification
```

### Horizon and extremal tiers

```@docs
kerr_geo_horizon_stable
kerr_geo_horizon_scatter
kerr_geo_extremal_family
kerr_geo_extremal
KerrGeoExtremalFamily
```

## Classification

```@docs
kerr_geo_classify
KerrGeoClassification
kerr_geo_components
KerrGeoRadialComponent
KerrGeoRadialEndpoint
kerr_geo_select_component
kerr_geo_classification_pipeline
kerr_geo_root_structure
```

## Classes and cases

```@docs
KERR_GEO_CLASSES
kerr_geo_class
kerr_geo_case_class
kerr_geo_case
KerrGeoCaseSpec
kerr_geo_case_catalog
kerr_geo_case_name
kerr_geo_case_symbol
kerr_geo_tier
HORIZON_TIER_IDS
EXTREMAL_TIER_IDS
kerr_geo_is_critical
kerr_geo_critical_role
CRITICAL_ROLES
```

## The Kerr metric and the potentials

```@docs
kerr_horizons
kerr_delta
kerr_metric_limit
kerr_energy_regime
kerr_energy_sign
kerr_radial_momentum
kerr_radial_potential
kerr_radial_coefficients
kerr_radial_polynomial
kerr_radial_derivatives
kerr_root_multiplicity_at
kerr_polar_z_potential
kerr_polar_theta_potential
kerr_polar_admissibility
kerr_polar_sector_candidates
kerr_axis_carter_q
kerr_axis_radial_potential
kerr_rstar
```

## APEX interface

```@docs
kerr_geo_constants_of_motion
kerr_geo_frequencies
kerr_geo_orbit
kerr_geo_stable
KerrGeoStable
kerr_geo_four_velocity
kerr_geo_radial_roots
kerr_geo_polar_roots
kerr_geo_orbit_type
kerr_geo_orbit_type_metadata
kerr_geo_isco
kerr_geo_isso
kerr_geo_ibso
kerr_geo_separatrix
```

## Finite-window interfaces

```@docs
kerr_geo_plunge
KerrGeoPlunge
radial_roots
polar_roots
classify_orbit
lambda_of_r
generic_plunge_orbit
generic_plunge_velocity
kerr_geo_capture
KerrGeoCapture
kerr_geo_scatter
KerrGeoScatter
```

## Diagnostics

```@docs
kerr_geo_diagnose
```
