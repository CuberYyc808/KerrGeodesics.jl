module KerrGeodesics

include("ConstantsOfMotion.jl")
using .ConstantsOfMotion
include("FourVelocity.jl")
using .FourVelocity
include("OrbitalFrequencies.jl")
using .OrbitalFrequencies
include("SpecialOrbits.jl")
using .SpecialOrbits
include("KerrGeoOrbit.jl")
using .KerrGeoOrbit

export Constants_Of_Motion,
        Four_Velocity,
        Orbital_Frequencies,
        Orbit_Type,
        Kerr_Geodesics,
        kerr_geo_emri

"""
    Constants_Of_Motion(a::Real, p::Real, e::Real, x::Real)
    returns Dict("E"=>_, "Lz"=>_, "Q"=>_) the three constants of motion for a particle on the geodesics of a Kerr black hole.
    "E" is energy per mass, "Lz" is angular momentum per mass, "Q" is Carter constant per mass square.
"""

Constants_Of_Motion(a::Real, p::Real, e::Real, x::Real) = kerr_geo_constants_of_motion(a, p, e, x)


"""
    Four_Velocity(a::Real, p::Real, e::Real, x::Real; initPhases=(0.0,0.0), Covariant=false, Parametrization="Mino")
    returns [ut(λ), ur(λ), uθ(λ), uϕ(λ)] as functions of Mino time λ.
    The index "initPhases" are the initial phases of the "r" and "θ" direction.
    The index "Covariant=false/true" determines whether contravariant or covariant form of the 4-velocity is returned.
"""

Four_Velocity(a::Real, p::Real, e::Real, x::Real; initPhases=(0.0, 0.0), Covariant=false, Parametrization="Mino") = kerr_geo_four_velocity(a, p, e, x; initPhases=initPhases, Covariant=Covariant, Parametrization=Parametrization)

"""
    Orbital_Frequencies(a::Real, p::Real, e::Real, x::Real; Time="Mino")
    returns Dict("ϒt", "ϒr", "ϒθ", "ϒϕ") if the parametrization is "Mino",
    or Dict("Ωr", "Ωθ", "Ωϕ") if the parametrization is "BoyerLindquist" or "Proper".
"""

Orbital_Frequencies(a::Real, p::Real, e::Real, x::Real; Time="Mino") = kerr_geo_frequencies(a, p, e, x; Time=Time)

"""
    Orbit_Type(a::Real, p::Real, e::Real, x::Real)
    returns the type of the orbit.
    For example ["Bound", "Eccentric", "Inclined"]
"""

Orbit_Type(a::Real, p::Real, e::Real, x::Real) = kerr_geo_orbit_type(a, p, e, x)

"""
    Kerr_Geodesics(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    returns a large association containing all the information of the set of parameter (a, p, e, x).
    One can access the trajectory through
        assco = KerrGeoOrbit(a, p, e, x; initPhases)
        t, r, θ, ϕ = assoc["Trajectory"]
    the BoyerLindquist coordinates are parametrized by Mino time "λ".
"""

Kerr_Geodesics(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0)) = kerr_geo_orbit(a, p, e, x; initPhases = initPhases)


struct KerrGeoEMRI
    OrbitalType::Vector{String}
    OrbitalParameters::NamedTuple
    ConstantsOfMotion::NamedTuple
    Parametrization::String
    Trajectory::NamedTuple
    InitialPhases::NamedTuple
    FourVelocity::NamedTuple
    Frequencies::NamedTuple
    CrossFunctions::NamedTuple
    DCrossFunctions::NamedTuple
end

function Base.show(io::IO, ::MIME"text/plain", kg::KerrGeoEMRI)
    println(io, "KerrGeoEMRI(")
    print(io, "    OrbitalType="); show(io, kg.OrbitalType); println(io, ",")
    print(io, "    OrbitalParameters="); show(io, kg.OrbitalParameters); println(io, ",")
    print(io, "    ConstantsOfMotion="); show(io, kg.ConstantsOfMotion); println(io, ",")
    print(io, "    Parametrization="); show(io, kg.Parametrization); println(io, ",")
    print(io, "    Trajectory="); show(io, kg.Trajectory); println(io, ",")
    print(io, "    InitialPhases="); show(io, kg.InitialPhases); println(io, ",")
    print(io, "    Frequencies="); show(io, kg.Frequencies); println(io, ",")
    print(io, ")")
end

function kerr_geo_emri(a::Real, p::Real, e::Real, x::Real; initPhases = (0.0, 0.0, 0.0, 0.0))
    # Orbital Type
    otype = Orbit_Type(a, p, e, x)

    # Constants of Motion
    com = Constants_Of_Motion(a, p, e, x)
    En = com["E"]
    L = com["Lz"]
    Q = com["Q"]

    # Trajectory
    KG = Kerr_Geodesics(a, p, e, x; initPhases = initPhases)
    t, r, θ, ϕ = KG["Trajectory"]
    # Frequencies
    freqs = Orbital_Frequencies(a, p, e, x; Time="Mino")
    ϒt = freqs["ϒt"]
    ϒr = freqs["ϒr"]
    ϒθ = freqs["ϒθ"]
    ϒϕ = freqs["ϒϕ"]
    # Cross functions
    Δtr = KG["CrossFunction"][1]
    Δtθ = KG["CrossFunction"][2]
    Δϕr = KG["CrossFunction"][3]
    Δϕθ = KG["CrossFunction"][4]
    # Derivatives of cross functions
    dtr = KG["DerivativesCrossFunction"][1]
    dtθ = KG["DerivativesCrossFunction"][2]
    dϕr = KG["DerivativesCrossFunction"][3]
    dϕθ = KG["DerivativesCrossFunction"][4]
    # Four-velocity
    velocity = Four_Velocity(a, p, e, x; initPhases=initPhases, Covariant=false, Parametrization="Mino")
    
    return KerrGeoEMRI(
        otype,
        (a=a, p=p, e=e, x=x),
        (E=En, Lz=L, Q=Q),
        "Mino",
        (t=t, r=r, θ=θ, ϕ=ϕ),
        (qr0=initPhases[1], qθ0=initPhases[2]),
        (ut=velocity[1], ur=velocity[2], uθ=velocity[3], uϕ=velocity[4]),
        (ϒt=ϒt, ϒr=ϒr, ϒθ=ϒθ, ϒϕ=ϒϕ),
        (Δtr=Δtr, Δtθ=Δtθ, Δϕr=Δϕr, Δϕθ=Δϕθ),
        (dtr=dtr, dtθ=dtθ, dϕr=dϕr, dϕθ=dϕθ)
    )
end



end