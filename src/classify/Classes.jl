# The six broad classes and the member tiers. Everything that names, orders or draws a class
# reads this table: case letter, Symbol, display name, KerrGeodesicFamily slot, catalogue order
# (the order of the rows) and plot colours (particle colour and four-stop trail palette).

"""
    KERR_GEO_CLASSES

The six orbit classes, in catalogue order. Each row is a NamedTuple with the class `symbol`
(`:stable`, `:critical`, `:plunge`, `:capture`, `:scatter`, `:trapped`), the case `letter`
(A, K, B, C, D, N), the display `name`, the `slot` of `KerrGeodesicFamily` that holds its
members, and the plot `color` and four-stop `palette` used by the example animations.
"""
const KERR_GEO_CLASSES = (
    (symbol=:stable, letter='A', name="Stable", slot=:Stable, color="#7ef0ff",
        palette=("#0b3d4f", "#1c8aa6", "#52d4e8", "#e8fdff")),
    (symbol=:critical, letter='K', name="Critical", slot=:Critical, color="#d9f36b",
        palette=("#2c3a0c", "#6f8f1c", "#c6e84a", "#f7ffd9")),
    (symbol=:plunge, letter='B', name="Plunge", slot=:Plunge, color="#ffc46b",
        palette=("#4a1606", "#b8420f", "#ff9a3c", "#fff1cf")),
    (symbol=:capture, letter='C', name="Capture", slot=:Capture, color="#ff7eb6",
        palette=("#44102e", "#b02a6c", "#ff6fae", "#ffe6f1")),
    (symbol=:scatter, letter='D', name="Scatter", slot=:Scatter, color="#b9a4ff",
        palette=("#1d1747", "#5a44c7", "#a591ff", "#f0ecff")),
    (symbol=:trapped, letter='N', name="Trapped", slot=:Trapped, color="#ff6b57",
        palette=("#4a0e0a", "#b3261e", "#ff6b57", "#ffe4df")),
)

"""
    kerr_geo_class(x)

Row of `KERR_GEO_CLASSES` for a class Symbol (`:critical`), a case letter (`'K'`) or a case ID
(`:K7`, `:D_H1`).
"""
function kerr_geo_class(x::Union{Symbol,Char})
    for row in KERR_GEO_CLASSES
        (x === row.symbol || x === row.letter) && return row
    end
    x isa Symbol && return kerr_geo_class(first(String(x)))
    throw(ArgumentError("$(repr(x)) is not a class, class letter or case ID."))
end

"""
    kerr_geo_case_class(id)

The class Symbol of a case or tier-member ID, read from its letter:
`kerr_geo_case_class(:K7) == :critical`, `kerr_geo_case_class(:D_H1) == :scatter`.
"""
function kerr_geo_case_class(id::Symbol)
    haskey(_CASE_BY_ID, id) || kerr_geo_tier(id) !== :primary ||
        throw(ArgumentError("$(id) is not a current case or member ID."))
    return kerr_geo_class(first(String(id))).symbol
end

# Members outside the primary numbering: H = horizon-root (P(r+) = 0, the horizon is a
# root of R), X = extremal Kerr (|a| = 1). Symbols use `_` (:D_H1), display names `-` ("D-H1").
"""
    HORIZON_TIER_IDS

`(:A_H1, :A_H2, :D_H1, :D_H2)`: the members of the horizon tier, which exist for
0 < |a| < 1 when P(r₊) = 0, so that the outer horizon is itself a root of R.
"""
const HORIZON_TIER_IDS = (:A_H1, :A_H2, :D_H1, :D_H2)

"""
    EXTREMAL_TIER_IDS

`(:A_X1, :A_X2, :B_X1, :B_X2, :C_X1, :C_X2, :C_X3, :C_X4, :D_X1, :D_X2)`: the members of
the extremal tier, which exist at |a| = 1 when P(r₊) = 2E − aLz = 0, so that the horizon is
a double or triple root of R.
"""
const EXTREMAL_TIER_IDS = (:A_X1, :A_X2, :B_X1, :B_X2, :C_X1, :C_X2, :C_X3, :C_X4,
    :D_X1, :D_X2)

"""
    kerr_geo_tier(id)

The tier of a member ID: `:horizon` for the H members, `:extremal` for the X members and
`:primary` for every case ID. A member built at |a| = 1 reports `Tier == :extremal` even
when it carries a primary case ID.
"""
kerr_geo_tier(id::Symbol) = id in HORIZON_TIER_IDS ? :horizon :
    id in EXTREMAL_TIER_IDS ? :extremal : :primary

"""
    kerr_geo_case_name(id)

The display name of a case or member ID: `kerr_geo_case_name(:D_H1) == "D-H1"`,
`kerr_geo_case_name(:K7) == "K7"`.
"""
kerr_geo_case_name(id::Symbol) =
    kerr_geo_tier(id) === :primary ? String(id) : replace(String(id), '_' => '-'; count=1)

"""
    kerr_geo_case_symbol(name)

The inverse of `kerr_geo_case_name`: `kerr_geo_case_symbol("D-H1") == :D_H1`.
"""
kerr_geo_case_symbol(name::Union{AbstractString,Symbol}) =
    Symbol(replace(String(name), '-' => '_'; count=1))
