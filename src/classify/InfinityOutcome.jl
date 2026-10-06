# Which closed-form E >= 1 construction applies to a set of constants.

function radial_roots_for_constants(a::Real, energy::Real, lz::Real, q::Real)
    coeffs = [-a^2 * q, 2 * (a * energy - lz)^2 + 2 * q,
              -(q + lz^2 - a^2 * _e2m1(energy)), 2, _e2m1(energy)]
    return _polynomial_roots(coeffs)
end

_real_roots(values; atol=_tol(real(float(eltype(values))), 1e-10)) =
    sort!([real(v) for v in values if abs(imag(v)) <= atol])

_root_class_symbol(n) = n == 4 ? :four_real : n == 3 ? :three_real :
    n == 2 ? :two_real_complex_pair : n == 1 ? :one_real : :other

"""
    _outcome_at_infinity(a, E, Lz, Q)

Classify constants for the finite-window scatter/capture constructors:
`formula` is `:parabolic_scatter` / `:hyperbolic_scatter` when the component
connected to infinity has an exterior turning point, `:parabolic_capture` /
`:hyperbolic_capture` when it reaches the future horizon, `:critical` otherwise, and
`:elliptic` for E < 1 (no motion to infinity). `energy_regime` is `:elliptic`, `:parabolic`
or `:hyperbolic`.
"""
function _outcome_at_infinity(a::Real, energy::Real, lz::Real, q::Real;
        atol=_tol(_float_type(a, energy, lz, q), 1e-10))
    rplus = _rplus(a)
    real_roots = _real_roots(radial_roots_for_constants(a, energy, lz, q); atol=atol)
    root_class = _root_class_symbol(length(real_roots))
    regime = kerr_energy_regime(energy)
    record(formula, outcome, reason) = (formula=formula, outcome=outcome,
        energy_regime=regime, root_class=root_class, roots=Tuple(real_roots), reason=reason)
    regime === :elliptic && return record(:elliptic, :plunge, "E < 1: no motion to infinity.")
    if any(r -> r > rplus + sqrt(atol), real_roots)
        return record(Symbol(regime, :_scatter), :scatter,
            "The component connected to infinity has an exterior turning point.")
    elseif kerr_radial_potential(a, energy, lz, q, rplus + 1e-7) >= -sqrt(atol)
        return record(Symbol(regime, :_capture), :capture,
            "No radial barrier between infinity and the future horizon.")
    end
    return record(:critical, :critical,
        "Neither a regular scatter nor a regular capture component.")
end
