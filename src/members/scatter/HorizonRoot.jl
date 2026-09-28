# Horizon-root scatter members (P(r+) = 0): D-H1 (E = 1), D-H2 (E > 1).

"""
    kerr_geo_horizon_scatter(a, E, Lz, Q; polar_sector=nothing, polar_phase=0.0,
                             polar_hemisphere=:north)

The horizon-root scatter member D-H1 (E = 1) or D-H2 (E > 1) for 0 < |a| < 1: P(r₊) = 0
makes r₊ a simple root of R, and the orbit comes in from infinity, turns at its single
exterior root (λ = 0, where t, φ, τ vanish) and returns to infinity.
"""
function kerr_geo_horizon_scatter(a::Real,energy::Real,lz::Real,q::Real;
        polar_sector=nothing,polar_phase::Real=0.0,
        polar_hemisphere::Symbol=:north)
    metric=kerr_metric_limit(a)
    metric in (:subextremal,:near_extremal) || throw(DomainError(
        a,"A subextremal horizon-root scatter orbit requires 0<|a|<1."))
    energy>=1 || throw(DomainError(
        energy,"A horizon-root scatter orbit requires E>=1."))
    horizons=kerr_horizons(a)
    pplus=kerr_radial_momentum(a,energy,lz,horizons.rplus)
    abs(pplus)<=2e-10*max(1.0,abs(energy),abs(lz)) || throw(DomainError(
        pplus,"A horizon-root scatter orbit requires P(r₊) = 0."))
    polar_admissibility=kerr_polar_admissibility(a,energy,lz,q)
    polar_admissibility.admissible || throw(DomainError(
        (energy,lz,q),"The polar potential is inadmissible."))
    structure=kerr_geo_root_structure(a,energy,lz,q)
    length(structure.horizon_coincident)==1 &&
        only(structure.horizon_coincident).multiplicity==1 || error(
        "A horizon-root scatter orbit requires a simple radial root on the outer horizon.")
    length(structure.exterior)==1 &&
        only(structure.exterior).multiplicity==1 || error(
        "The horizon-root scatter component requires one simple outer turn.")
    selected=polar_sector===nothing ? _constants_polar_sector(energy,lz,q) : polar_sector
    polar=_polar_solution(a,energy,lz,q,selected,float(polar_phase);
        hemisphere=polar_hemisphere)
    parabolic=energy==1
    return _scatter_member(a,energy,lz,q,parabolic ? :D_H1 : :D_H2,polar,
        Tuple(root.radius for root in structure.real_roots);
        formula_family=parabolic ? :HC_FF07 : :HC_FF12)
end
