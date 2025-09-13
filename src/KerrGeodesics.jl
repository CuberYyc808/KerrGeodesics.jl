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
        Kerr_Geodesics

# Integration of some important functions.

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

end
