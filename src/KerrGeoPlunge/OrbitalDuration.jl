module OrbitalDuration

include("OrbitClass.jl")
using .OrbitClass
using Elliptic

export lambda_of_r

function lambda_of_r_real(a, E, roots)
    r4, r3, r2, r1 = roots
    function λ_of_r(r)
        if r4 <= r <= r3
            yr = sqrt((r - r3) / (r - r2) * (r2 - r4) / (r3 - r4))
            kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
            return 2 * Elliptic.F(asin(yr), kr) / sqrt((1 - E^2) * (r1 - r3) * (r2 - r4))
        else 
            @info("r = $r is out of the deeply bound plunge region between r4 = $r4 and r3 = $r3.")
        end
    end
    Λr_max = λ_of_r(r4)
    rp = 1 + sqrt(1 - a^2)
    Λr_tilde_max = λ_of_r(rp)
    return Λr_max, Λr_tilde_max, λ_of_r
end

function lambda_of_r_complex(a, E, roots)
    r1, r2, A, B = roots
    function λ_of_r(r)
        if r2 <= r <= r1
            yr = (B * (r1 - r) - A * (r - r2)) / (B * (r1 - r) + A * (r - r2))
            kr = ((r1 - r2)^2 - (A - B)^2) / (4 * A * B)
            return Elliptic.F(π/2 + asin(yr), kr) / sqrt((1 - E^2) * A * B)
        else 
            @info("r = $r is out of the plunge region between r2 = $r2 and r1 = $r1.")
        end
    end
    Λr_max = λ_of_r(r2)
    rp = 1 + sqrt(1 - a^2)
    Λr_tilde_max = λ_of_r(rp)
    return Λr_max, Λr_tilde_max, λ_of_r
end

function lambda_of_r_real2(a, E, roots)
    r4, r3, r2, r1 = roots
    rp = 1 + sqrt(1 - a^2)
    kr = (r1 - r2) / (r1 - r3) * (r3 - r4) / (r2 - r4)
    ξr = sqrt((1 - E^2) * (r1 - r3) * (r2 - r4)) / 2
    k_complete = Elliptic.K(kr)
    function λ_of_r(r)
        if rp <= r <= r1
            u = (r - r2) * (r1 - r3) / ((r1 - r2) * (r - r3))
            u = clamp(u, 0.0, 1.0)
            return (k_complete - Elliptic.F(asin(sqrt(u)), kr)) / ξr
        else
            @info("r = $r is out of the Real2 exterior plunge region between r+ = $rp and r1 = $r1.")
        end
    end
    Λr_max = λ_of_r(rp)
    Λr_tilde_max = Λr_max
    return Λr_max, Λr_tilde_max, λ_of_r
end

"""
    lambda_of_r(a, E, L, Q)

Return branch-specific support radii and a radius-to-Mino-time map used by the
plunge wrapper.
"""
function lambda_of_r(a, E, L, Q)
    roots, cf = classify_orbit(a, E, L, Q)
    if cf == "Real1"
        return lambda_of_r_real(a, E, roots)
    elseif cf == "Complex"
        return lambda_of_r_complex(a, E, roots)
    elseif cf == "Real2"
        return lambda_of_r_real2(a, E, roots)
    else 
        @info("The orbit is classified as $cf which has not been included in the implementation.")
    end
end

end
