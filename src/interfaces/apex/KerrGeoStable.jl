# APEX reference API: `kerr_geo_stable`, the KerrGeoStable record of a Stable orbit from (a, p, e, x).

"""
    kerr_geo_stable(a, p, e, x; initPhases=(0.0, 0.0, 0.0, 0.0))

The stable orbit with APEX parameters `(a, p, e, x)` as a `KerrGeoStable` record. The
initial phases `(qt0, qr0, qθ0, qφ0)` shift t, the radial phase, the polar phase and φ at
λ = 0; with zero phases, λ = 0 is at periapsis and at the northern polar turning point.
"""
function kerr_geo_stable(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    # Orbital Type
    otype = kerr_geo_orbit_type(a, p, e, x)

    # Constants of Motion
    com = kerr_geo_constants_of_motion(a, p, e, x)
    En = com["E"]
    L = com["Lz"]
    Q = com["Q"]

    # Trajectory
    KG = kerr_geo_orbit(a, p, e, x; initPhases = initPhases)
    t, r, θ, ϕ = KG["Trajectory"]
    # Frequencies
    freqs = kerr_geo_frequencies(a, p, e, x; Time="Mino")
    ϒt = freqs["ϒt"]
    ϒr = freqs["ϒr"]
    ϒθ = freqs["ϒθ"]
    ϒϕ = freqs["ϒϕ"]
    # Cross functions
    if KG["CrossFunction"] !== nothing
        Δtr = KG["CrossFunction"][1]
        Δtθ = KG["CrossFunction"][2]
        Δϕr = KG["CrossFunction"][3]
        Δϕθ = KG["CrossFunction"][4]
    else
        Δtr = nothing
        Δtθ = nothing
        Δϕr = nothing
        Δϕθ = nothing
    end
    # Derivatives of cross functions
    if KG["DerivativesCrossFunction"] !== nothing
        dtr = KG["DerivativesCrossFunction"][1]
        dtθ = KG["DerivativesCrossFunction"][2]
        dϕr = KG["DerivativesCrossFunction"][3]
        dϕθ = KG["DerivativesCrossFunction"][4]
    else
        dtr = nothing
        dtθ = nothing
        dϕr = nothing
        dϕθ = nothing
    end
    # Four-velocity
    ut, ur, uθ, uϕ = KG["FourVelocity"]
    
    return KerrGeoStable(
        otype,
        (a=a, p=p, e=e, x=x),
        (E=En, Lz=L, Q=Q),
        "Mino",
        (t=t, r=r, θ=θ, ϕ=ϕ),
        (qt0 = initPhases[1], qr0=initPhases[2], qθ0=initPhases[3], qϕ0=initPhases[4]),
        (ut=ut, ur=ur, uθ=uθ, uϕ=uϕ),
        (ϒt=ϒt, ϒr=ϒr, ϒθ=ϒθ, ϒϕ=ϒϕ),
        (Δtr=Δtr, Δtθ=Δtθ, Δϕr=Δϕr, Δϕθ=Δϕθ),
        (dtr=dtr, dtθ=dtθ, dϕr=dϕr, dϕθ=dϕθ)
    )
end
