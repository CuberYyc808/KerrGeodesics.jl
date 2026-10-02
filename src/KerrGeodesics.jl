# Layout (by responsibility): core/ (Chebyshev spectral tools, metric functions, the polar and
# radial coordinate engines) -> classify/ (class table, case table and classifiers) -> models/
# (radial and polar motion models) -> members/ (the member type and its shared assembly, one
# directory per class: Stable, Critical, Plunge, Capture, Scatter, Trapped, plus the
# exact-extremal members) -> interfaces/ (the APEX-parameter API: constants of motion,
# frequencies, `kerr_geo_orbit`, `kerr_geo_plunge`) -> family/ (`kerr_geodesic`) -> Diagnostics.

"""
    KerrGeodesics

Timelike geodesics of the Kerr spacetime in Mino time. `kerr_geodesic(a, (E, Lz, Q))` returns
every orbit the constants of motion allow outside the black hole, each with its coordinates,
four-velocity and, for stable orbits, frequencies as functions of Mino time λ; units are
G = c = M = 1.
"""
module KerrGeodesics

using Elliptic
using Polynomials
using QuadGK
using Roots

include("core/Spectral.jl")
include("core/Elliptic.jl")
include("core/Series.jl")
include("core/ExactInputArithmetic.jl")
include("core/Metric.jl")
include("core/Degeneracy.jl")
include("core/PolarEngine.jl")
include("core/RadialEngine.jl")
include("classify/Classes.jl")
include("classify/Cases.jl")
include("classify/Classifier.jl")
include("classify/InfinityOutcome.jl")
include("models/RadialPrimitives.jl")
include("models/InteriorRepeatedRoots.jl")
include("models/SimpleRootModels.jl")
include("models/HorizonRelative.jl")
include("models/FourRealRoots.jl")
include("models/ThreeRealRoots.jl")
include("models/WindowPolarMotion.jl")
include("models/PolarSolutions.jl")
include("models/ApexConstants.jl")
include("members/Kinematics.jl")
include("members/Component.jl")
include("members/Assembly.jl")
include("members/critical/Models.jl")
include("members/critical/CriticalMembers.jl")
include("members/stable/Landmarks.jl")
include("members/stable/StableCases.jl")
include("members/stable/StableOrbit.jl")
include("members/stable/StableComponent.jl")
include("members/stable/HorizonRoot.jl")
include("members/plunge/AxisInfall.jl")
include("members/plunge/PlungeCases.jl")
include("members/scatter/FiniteWindow.jl")
include("members/scatter/ScatterCases.jl")
include("members/scatter/HorizonRoot.jl")
include("members/capture/FiniteWindow.jl")
include("members/capture/FourComplex.jl")
include("members/capture/CaptureCases.jl")
include("members/trapped/TrappedCases.jl")
include("members/extremal/Extremal.jl")
include("members/extremal/EngineMembers.jl")
include("interfaces/apex/ConstantsOfMotion.jl")
include("interfaces/apex/OrbitalFrequencies.jl")
include("interfaces/apex/FourVelocity.jl")
include("interfaces/apex/KerrGeoOrbit.jl")
include("interfaces/apex/KerrGeoStable.jl")
include("interfaces/plunge/OrbitClass.jl")
include("interfaces/plunge/OrbitalDuration.jl")
include("interfaces/plunge/FourVelocity.jl")
include("interfaces/plunge/PlungeOrbit.jl")
include("interfaces/plunge/HorizonSeries.jl")
include("interfaces/plunge/InitialConditions.jl")
include("interfaces/plunge/PlungeApi.jl")
include("family/Family.jl")
include("Diagnostics.jl")
include("Precompile.jl")

export kerr_geo_constants_of_motion,
        kerr_delta,
        kerr_horizons,
        kerr_metric_limit,
        kerr_energy_regime,
        kerr_energy_sign,
        kerr_radial_momentum,
        kerr_axis_carter_q,
        kerr_axis_radial_potential,
        kerr_radial_coefficients,
        kerr_radial_polynomial,
        kerr_radial_potential,
        kerr_radial_derivatives,
        kerr_root_multiplicity_at,
        kerr_polar_theta_potential,
        kerr_polar_z_potential,
        kerr_polar_admissibility,
        kerr_polar_sector_candidates,
        KerrGeoCaseSpec,
        KerrGeoRadialEndpoint,
        KerrGeoRadialComponent,
        KerrGeoClassification,
        kerr_geo_case,
        kerr_geo_case_catalog,
        kerr_geo_case_name,
        kerr_geo_case_symbol,
        kerr_geo_case_class,
        KERR_GEO_CLASSES,
        kerr_geo_class,
        kerr_geo_tier,
        HORIZON_TIER_IDS,
        EXTREMAL_TIER_IDS,
        kerr_geo_is_critical,
        kerr_geo_critical_component,
        kerr_geo_members,
        kerr_geo_critical_role,
        CRITICAL_ROLES,
        kerr_geo_root_structure,
        kerr_geo_components,
        kerr_geo_classify,
        kerr_geo_select_component,
        kerr_geo_classification_pipeline,
        kerr_geo_stability_metadata,
        KerrGeoComponent,
        KerrGeoCriticalComponent,
        kerr_geo_member_class,
        kerr_geo_sample,
        kerr_geo_four_velocity,
        kerr_geo_frequencies,
        kerr_geo_radial_roots,
        kerr_geo_polar_roots,
        kerr_geo_orbit_type,
        kerr_geo_orbit_type_metadata,
        kerr_geo_separatrix,
        kerr_geo_isco,
        kerr_geo_ibso,
        kerr_geo_isso,
        kerr_geo_orbit,
        kerr_geo_stable,
        KerrGeoStable,
        KerrGeoStableComponent,
        kerr_geo_stable_component,
        KerrGeoExtremalFamily,
        kerr_geo_extremal,
        kerr_geo_extremal_family

export radial_roots, polar_roots, classify_orbit, 
        lambda_of_r, 
        generic_plunge_velocity,
        generic_plunge_orbit,
        kerr_rstar,
        KerrGeoPlunge,
        KerrGeodesicFamily,
        kerr_geo_plunge,
        KerrGeoPlungeComponent,
        kerr_geo_plunge_component,
        kerr_geo_critical_spherical,
        kerr_geo_plunge_axis_infall,
        kerr_geo_critical_plunge,
        kerr_geo_capture_axis_infall,
        kerr_geo_capture_vortical,
        kerr_geo_capture_constant_latitude,
        kerr_geo_capture_equator_attractive,
        kerr_geo_capture_four_complex,
        kerr_geo_critical_homoclinic,
        KerrGeoScatter,
        KerrGeoScatterComponent,
        kerr_geo_scatter_component,
        KerrGeoCapture,
        KerrGeoCaptureComponent,
        kerr_geo_capture_component,
        KerrGeoTrappedComponent,
        KerrGeoTrappedClassification,
        kerr_geo_trapped,
        kerr_geo_trapped_case,
        kerr_geo_trapped_classify,
        kerr_geo_scatter_asymptotic_diagnostics,
        kerr_geo_scatter_asymptotic_state,
        kerr_geo_scatter,
        kerr_geo_capture,
        kerr_geodesic

export kerr_geo_horizon_stable,
       kerr_geo_horizon_scatter

export kerr_geo_diagnose

end
