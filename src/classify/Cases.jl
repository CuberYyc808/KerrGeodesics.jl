# Case table: one KerrGeoCaseSpec per primary case (A1-A2, K1-K11, B1-B9, C1-C12, D1-D2,
# N1-N6), and the records the classifier returns (endpoints, components, classification).

"""
    KerrGeoCaseSpec

The definition of one primary case, as returned by `kerr_geo_case(id)`: its `CaseId`, the
`BroadClass`, the `EnergyRegime` (the sign of E² − 1) and `EnergySign`, the `Degree` and
`LeadingSign` of R, the root structure (`RootTopology`, the ordering of the roots and the
horizon in `RootOrdering`, their `RootMultiplicities`, and `RootsAboveHorizon`), the
`AllowedInterval` of r, the `PastEndpoint` and `FutureEndpoint` of the motion, its
`RadialOrientation`, the cases that share its constants (`PairedCaseIds`,
`FamilyMemberCaseIds`), and the `FormulaFamily` of its closed-form r(λ).
"""
struct KerrGeoCaseSpec
    CaseId::Symbol
    BroadClass::Symbol
    EnergyRegime::Symbol
    EnergySign::Int
    Degree::Int
    LeadingSign::Int
    RootTopology::Symbol
    RootOrdering::String
    RootMultiplicities::Tuple{Vararg{Int}}
    RootsAboveHorizon::Int
    AllowedInterval::String
    PastEndpoint::Symbol
    FutureEndpoint::Symbol
    RadialOrientation::Symbol
    PairedCaseIds::Tuple{Vararg{Symbol}}
    FamilyMemberCaseIds::Tuple{Vararg{Symbol}}
    FormulaFamily::Symbol
end

"""
    KerrGeoRadialEndpoint

One end of an allowed radial interval: its `Kind` (`:radial_root`, `:outer_horizon` or
`:infinity`), `Radius`, whether the interval includes it (`Included`: true for a turning
point; false for the horizon, infinity and a repeated root approached asymptotically), and
the `Multiplicity` of the root there. `T` is the floating-point type of the constants.
"""
struct KerrGeoRadialEndpoint{T<:Real}
    Kind::Symbol
    Radius::T
    Included::Bool
    Multiplicity::Int
end

"""
    KerrGeoRadialComponent

One region of radial motion found by `kerr_geo_classify`: an interval between two
`KerrGeoRadialEndpoint`s (`LowerEndpoint`, `UpperEndpoint`) on which R ≥ 0, or a repeated
root the orbit can sit on. It carries the `CaseId` and `BroadClass` assigned to it, the
`EnergyRegime`, its `Connectivity` and `RadialOrientation`, the `FormulaFamily`, the cases
sharing its constants (`PairedCaseIds`, `FamilyMemberCaseIds`), the `PolarSector` of the
accompanying polar motion, and `Tags`, `SupportStatus` and `Metadata` (roots, radial
derivatives at the endpoints).
"""
struct KerrGeoRadialComponent{T<:Real}
    CaseId::Union{Nothing,Symbol}
    BroadClass::Symbol
    EnergyRegime::Symbol
    LowerEndpoint::KerrGeoRadialEndpoint{T}
    UpperEndpoint::KerrGeoRadialEndpoint{T}
    Connectivity::Symbol
    RadialOrientation::Symbol
    FormulaFamily::Symbol
    PairedCaseIds::Tuple{Vararg{Symbol}}
    FamilyMemberCaseIds::Tuple{Vararg{Symbol}}
    PolarSector::Symbol
    Tags::Tuple{Vararg{Symbol}}
    SupportStatus::Symbol
    Metadata::NamedTuple
end

"""
    KerrGeoClassification

The result of `kerr_geo_classify`: the input `Parameters` and `ConstantsOfMotion`, the
`EnergyRegime` and `MetricLimit`, the complex `Roots` of R, every admitted
`KerrGeoRadialComponent` (`Components`) with their `CaseIds`, the cases the constants rule
out (`ExcludedCaseIds`), the `PolarMetadata`, `Tags`, the case chosen by the selection
keywords (`SelectedCase`, `SelectionHint`) and a `Status`.
"""
struct KerrGeoClassification{T<:Real}
    Parameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    EnergyRegime::Symbol
    MetricLimit::Symbol
    Roots::Vector{Complex{T}}
    Components::Vector{KerrGeoRadialComponent{T}}
    CaseIds::Tuple{Vararg{Symbol}}
    ExcludedCaseIds::Tuple{Vararg{Symbol}}
    PolarMetadata::NamedTuple
    Tags::Tuple{Vararg{Symbol}}
    SelectedCase::Union{Nothing,Symbol}
    SelectionHint::NamedTuple
    Status::NamedTuple
end
KerrGeoClassification(parameters, constants, regime, limit, roots::Vector{Complex{T}},
    components::AbstractVector, rest...) where {T} =
    KerrGeoClassification{T}(parameters, constants, regime, limit, roots, components, rest...)

# catalogue order of a case ID: class order of KERR_GEO_CLASSES, then number
_case_order(id::Symbol) = (findfirst(row -> row.letter == first(String(id)), KERR_GEO_CLASSES),
    parse(Int, String(id)[2:end]))

# BroadClass comes from the case letter (Classes.jl); FamilyMemberCaseIds are all cases with the
# same constants, the case itself included.
function _spec(id, energy, degree, leading, topology, ordering,
        multiplicities, above, interval, past, future, orientation, pairs, family)
    family_members = isempty(pairs) ? () : Tuple(sort([id, pairs...]; by=_case_order))
    class = kerr_geo_class(first(String(id))).symbol
    return KerrGeoCaseSpec(
        id, class, energy, class === :trapped ? -1 : 1, degree, leading, topology, ordering,
        Tuple(multiplicities), above, interval, past, future, orientation,
        Tuple(pairs), family_members, family,
    )
end

const _CASE_SPECS = KerrGeoCaseSpec[
    _spec(:A1, :elliptic, 4, -1, :four_real_simple,
        "r1<rplus<r2<r3<r4", (1,1,1,1), 3, "[r3,r4]",
        :finite_turning_point, :finite_turning_point, :libration, (:B1,), :FF01),
    _spec(:A2, :elliptic, 4, -1, :four_real_outer_double,
        "r1<rplus<r2<r3=r4", (1,1,2), 3, "{r3}",
        :stable_repeated_root, :stable_repeated_root, :constant_radius, (:B2,), :FF02),
    _spec(:K1, :elliptic, 4, -1, :four_real_outer_triple,
        "r1<rplus<r2=r3=r4", (1,3), 3, "{r2}",
        :marginal_repeated_root, :marginal_repeated_root, :constant_radius, (:K2,), :FF03),
    _spec(:K2, :elliptic, 4, -1, :four_real_outer_triple,
        "r1<rplus<r2=r3=r4", (1,3), 3, "(rplus,r2)",
        :marginal_repeated_root, :future_horizon, :inward, (:K1,), :FF03),
    _spec(:K3, :elliptic, 4, -1, :four_real_unstable_double,
        "r1<rplus<rc=rc<ra", (1,2,1), 3, "{rc}",
        :unstable_repeated_root, :unstable_repeated_root, :constant_radius, (:K4, :K5), :FF17),
    _spec(:K4, :elliptic, 4, -1, :four_real_unstable_double,
        "r1<rplus<rc=rc<ra", (1,2,1), 3, "(rc,ra]",
        :unstable_repeated_root, :unstable_repeated_root, :homoclinic, (:K3, :K5), :FF17),
    _spec(:K5, :elliptic, 4, -1, :four_real_unstable_double,
        "r1<rplus<rc=rc<ra", (1,2,1), 3, "(rplus,rc)",
        :unstable_repeated_root, :future_horizon, :inward, (:K3, :K4), :FF17),
    _spec(:K6, :parabolic, 3, 1, :three_real_outer_double,
        "r1<rplus<r2=r3", (1,2), 2, "{r2}",
        :unstable_repeated_root, :unstable_repeated_root, :constant_radius, (:K7, :K8), :FF06),
    _spec(:K7, :parabolic, 3, 1, :three_real_outer_double,
        "r1<rplus<r2=r3", (1,2), 2, "(r3,inf)",
        :past_infinity, :unstable_repeated_root, :inward_asymptotic, (:K6, :K8), :FF06),
    _spec(:K8, :parabolic, 3, 1, :three_real_outer_double,
        "r1<rplus<r2=r3", (1,2), 2, "(rplus,r2)",
        :unstable_repeated_root, :future_horizon, :inward, (:K6, :K7), :FF06),
    _spec(:K9, :hyperbolic, 4, 1, :four_real_outer_double,
        "r1<r2<rplus<r3=r4", (1,1,2), 2, "{r3}",
        :unstable_repeated_root, :unstable_repeated_root, :constant_radius, (:K10, :K11),
        :FF11),
    _spec(:K10, :hyperbolic, 4, 1, :four_real_outer_double,
        "r1<r2<rplus<r3=r4", (1,1,2), 2, "(r4,inf)",
        :past_infinity, :unstable_repeated_root, :inward_asymptotic, (:K9, :K11), :FF11),
    _spec(:K11, :hyperbolic, 4, 1, :four_real_outer_double,
        "r1<r2<rplus<r3=r4", (1,1,2), 2, "(rplus,r3)",
        :unstable_repeated_root, :future_horizon, :inward, (:K9, :K10), :FF11),
    _spec(:B1, :elliptic, 4, -1, :four_real_simple,
        "r1<rplus<r2<r3<r4", (1,1,1,1), 3, "(rplus,r2]",
        :finite_turning_point, :future_horizon, :inward, (:A1,), :FF01),
    _spec(:B2, :elliptic, 4, -1, :four_real_outer_double,
        "r1<rplus<r2<r3=r4", (1,1,2), 3, "(rplus,r2]",
        :finite_turning_point, :future_horizon, :inward, (:A2,), :FF02),
    _spec(:B3, :elliptic, 4, -1, :four_real_single_exterior,
        "r1<r2<r3<rplus<r4", (1,1,1,1), 1, "(rplus,r4]",
        :finite_turning_point, :future_horizon, :inward, (), :FF04),
    _spec(:B4, :elliptic, 4, -1, :two_real_complex_pair,
        "r1<rplus<r2;complex_pair", (1,1), 1, "(rplus,r2]",
        :finite_turning_point, :future_horizon, :inward, (), :FF05),
    _spec(:B5, :parabolic, 3, 1, :three_real_simple,
        "r1<rplus<r2<r3", (1,1,1), 2, "(rplus,r2]",
        :finite_turning_point, :future_horizon, :inward, (:D1,), :FF07),
    _spec(:B6, :hyperbolic, 4, 1, :four_real_simple,
        "r1<r2<rplus<r3<r4", (1,1,1,1), 2, "(rplus,r3]",
        :finite_turning_point, :future_horizon, :inward, (:D2,), :FF12),
    _spec(:B7, :elliptic, 4, -1, :interior_double_then_simple,
        "rd=rd<s<rplus<ro", (2,1,1), 1, "(rplus,ro]",
        :finite_turning_point, :future_horizon, :inward, (), :ER01),
    _spec(:B8, :elliptic, 4, -1, :interior_simple_then_double,
        "s<rd=rd<rplus<ro", (1,2,1), 1, "(rplus,ro]",
        :finite_turning_point, :future_horizon, :inward, (), :ER02),
    _spec(:B9, :elliptic, 4, -1, :interior_triple,
        "rd=rd=rd<rplus<ro", (3,1), 1, "(rplus,ro]",
        :finite_turning_point, :future_horizon, :inward, (), :ER03),
    _spec(:C1, :parabolic, 3, 1, :one_real_complex_pair,
        "r1<rplus;complex_pair", (1,), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :FF08),
    _spec(:C2, :parabolic, 3, 1, :three_real_all_below_horizon,
        "r1<r2<r3<rplus", (1,1,1), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :FF09),
    _spec(:C3, :hyperbolic, 4, 1, :two_real_complex_pair_below_horizon,
        "r1<r2<rplus;complex_pair", (1,1), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :FF13),
    _spec(:C4, :hyperbolic, 4, 1, :four_real_all_below_horizon,
        "r1<r2<r3<r4<rplus", (1,1,1,1), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :FF14),
    _spec(:C5, :hyperbolic, 4, 1, :four_complex_no_real_roots,
        "two_complex_conjugate_pairs;no_real_roots", (), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :FF18),
    _spec(:C6, :parabolic, 3, 1, :interior_double_then_simple,
        "rd=rd<s<rplus", (2,1), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER04),
    _spec(:C7, :parabolic, 3, 1, :interior_simple_then_double,
        "s<rd=rd<rplus", (1,2), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER05),
    _spec(:C8, :parabolic, 3, 1, :interior_triple,
        "rd=rd=rd<rplus", (3,), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER06),
    _spec(:C9, :hyperbolic, 4, 1, :interior_simple_double_simple,
        "s1<rd=rd<s2<rplus", (1,2,1), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER07),
    _spec(:C10, :hyperbolic, 4, 1, :interior_two_simple_then_double,
        "s1<s2<rd=rd<rplus", (1,1,2), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER08),
    _spec(:C11, :hyperbolic, 4, 1, :interior_double_plus_complex_pair,
        "rd=rd<rplus;complex_pair", (2,), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER09),
    _spec(:C12, :hyperbolic, 4, 1, :interior_simple_then_triple,
        "s<rd=rd=rd<rplus", (1,3), 0, "(rplus,inf)",
        :past_infinity, :future_horizon, :inward, (), :ER10),
    _spec(:D1, :parabolic, 3, 1, :three_real_simple,
        "r1<rplus<r2<r3", (1,1,1), 2, "[r3,inf)",
        :past_infinity, :future_infinity, :inbound_turn_outbound, (:B5,), :FF07),
    _spec(:D2, :hyperbolic, 4, 1, :four_real_simple,
        "r1<r2<rplus<r3<r4", (1,1,1,1), 2, "[r4,inf)",
        :past_infinity, :future_infinity, :inbound_turn_outbound, (:B6,), :FF12),
    # Class N (E < 0): past horizon -> turning point -> future horizon inside the ergoregion;
    # one case per radial root structure: case Nk is radial disposition NFD0k
    # (`Status.disposition_id`)
    _spec(:N1, :elliptic, 4, -1, :four_real_simple,
        "r1<rplus<r2<r3<r4", (1,1,1,1), 3, "(rplus,r2]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF01),
    _spec(:N2, :elliptic, 4, -1, :four_real_outer_double,
        "r1<rplus<r2<r3=r4", (1,1,2), 3, "(rplus,r2]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF02),
    _spec(:N3, :elliptic, 4, -1, :four_real_single_exterior,
        "r1<r2<r3<rplus<r4", (1,1,1,1), 1, "(rplus,r4]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF04),
    _spec(:N4, :elliptic, 4, -1, :two_real_complex_pair,
        "r1<rplus<r2;complex_pair", (1,1), 1, "(rplus,r2]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF05),
    _spec(:N5, :parabolic, 3, 1, :three_real_simple,
        "r1<rplus<r2<r3", (1,1,1), 2, "(rplus,r2]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF07),
    _spec(:N6, :hyperbolic, 4, 1, :four_real_simple,
        "r1<r2<rplus<r3<r4", (1,1,1,1), 2, "(rplus,r3]",
        :past_horizon, :future_horizon, :outbound_turn_inbound, (), :FF12),
]

const _CASE_BY_ID = Dict(spec.CaseId => spec for spec in _CASE_SPECS)

# the six radial dispositions of Class N, one per case (NFD0k = Nk)
const _DISPOSITION_CASE_IDS = Dict(Symbol("NFD0", k) => Symbol("N", k) for k in 1:6)
const TRAPPED_CASE_IDS = (:N1, :N2, :N3, :N4, :N5, :N6)
# Class N case of a radial disposition (NFD0k -> Nk)
_trapped_case_id(disposition::Symbol) = _DISPOSITION_CASE_IDS[disposition]

"""
    kerr_geo_case(id) -> KerrGeoCaseSpec

Specification of a primary case (`:A1`, `:K7`, `:N3`, ...): broad class, energy regime and
sign, degree and leading sign of R, root topology, ordering and multiplicities, allowed radial
interval, past and future endpoints, radial orientation, paired and same-constants cases, and
formula family. H- and X-tier members have no entry.
"""
function kerr_geo_case(case_id::Symbol)
    haskey(_CASE_BY_ID, case_id) ||
        error("Unknown case $(case_id). kerr_geo_case_catalog() lists the cases.")
    return _CASE_BY_ID[case_id]
end

"""
    kerr_geo_case_catalog()

The table of primary cases, a `Dict` from case ID to `KerrGeoCaseSpec` (a copy; changing it
does not affect the package).
"""
kerr_geo_case_catalog() = Dict(key => value for (key, value) in _CASE_BY_ID)
