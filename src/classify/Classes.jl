# The six broad classes and the member tiers. Everything that names, orders or draws a class
# reads this table: case letter, Symbol, display name, KerrGeodesicFamily slot, catalogue order
# (the order of the rows) and plot colours (particle colour and four-stop trail palette).

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

"""Broad class Symbol of a case or tier-member ID, from its letter: `kerr_geo_case_class(:K7) == :critical`."""
function kerr_geo_case_class(id::Symbol)
    haskey(_CASE_BY_ID, id) || kerr_geo_tier(id) !== :primary ||
        throw(ArgumentError("$(id) is not a current case or member ID."))
    return kerr_geo_class(first(String(id))).symbol
end

# Members outside the primary numbering: H = horizon-root (P(r+) = 0, the horizon is a
# root of R), X = extremal Kerr (|a| = 1). Symbols use `_` (:D_H1), display names `-` ("D-H1").
const HORIZON_TIER_IDS = (:A_H1, :A_H2, :D_H1, :D_H2)
const EXTREMAL_TIER_IDS = (:A_X1, :A_X2, :B_X1, :B_X2, :C_X1, :C_X2, :C_X3, :C_X4,
    :D_X1, :D_X2)

"""Tier of a member ID: `:horizon` (H), `:extremal` (X) or `:primary`."""
kerr_geo_tier(id::Symbol) = id in HORIZON_TIER_IDS ? :horizon :
    id in EXTREMAL_TIER_IDS ? :extremal : :primary

"""Display name of a case or member ID: `:D_H1 -> "D-H1"`, `:K7 -> "K7"`."""
kerr_geo_case_name(id::Symbol) =
    kerr_geo_tier(id) === :primary ? String(id) : replace(String(id), '_' => '-'; count=1)

"""Inverse of `kerr_geo_case_name`: `"D-H1" -> :D_H1`."""
kerr_geo_case_symbol(name::Union{AbstractString,Symbol}) =
    Symbol(replace(String(name), '-' => '_'; count=1))
