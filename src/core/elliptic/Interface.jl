# Elliptic integrals and Jacobi functions in the parameter convention m of the APEX and plunge
# reference interfaces. NaN propagates; a parameter outside the domain of a form raises a
# DomainError: m > 1 for K and E (K(1) = ∞, E(1) = 1), m ∉ [0, 1] for Π and the incomplete
# integrals and the amplitude. sn, cn, dn take any real m.

function _check_parameter(m)
    (m < 0 || m > 1) && throw(DomainError(m, "argument m not in [0,1]"))
    return m
end
function _check_complete_parameter(m)
    m > 1 && throw(DomainError(m, "argument m not <= 1"))
    return m
end

# K(m), E(m), D(m) = (K − E)/m, Π(n|m)
function _K(m)
    isnan(m) && return float(m)
    _check_complete_parameter(m) == 1 && return float(oftype(m, Inf))
    return _ellip_k(1 - m)
end
function _E(m)
    isnan(m) && return float(m)
    _check_complete_parameter(m) == 1 && return one(float(m))
    return _ellip_e_complete(1 - m)
end
_D(m) = _ellip_d_complete(1 - m)
_Pi(n, m) = isnan(n) || isnan(m) ? float(n + m) : _complete_pi(n, 1 - n, _check_parameter(m))
_complete_pi(n, n1, m) = _ellip_pi_complete(1 - m, n, n1)

# F(φ|m), E(φ|m), D(φ|m), Π(n; φ|m) for any real amplitude
_F(φ, m) = isnan(φ) || isnan(m) ? float(φ + m) : _ellip_f(φ, 1 - _check_parameter(m))
function _E(φ, m)
    (isnan(φ) || isnan(m)) && return float(φ + m)
    m1 = 1 - _check_parameter(m)
    return _legendre((s, c) -> _ellip_e(s, c, m1), () -> _E(m), φ)
end
_D(φ, m) = _ellip_d(φ, 1 - m)
_Pi(n, φ, m) = isnan(n) || isnan(φ) || isnan(m) ? float(n + φ + m) :
    _ellip_pi(φ, 1 - _check_parameter(m), n, 1 - n)

# Jacobi functions: one parameter record per orbit, sn, cn, dn for any real m, the amplitude
# for 0 ≤ m ≤ 1
_sn(u, J::_JacobiParameter) = _jacobi_sncndn(u, J)[1]
_cn(u, J::_JacobiParameter) = _jacobi_sncndn(u, J)[2]
_dn(u, J::_JacobiParameter) = _jacobi_sncndn(u, J)[3]
function _am(u, J::_JacobiParameter)
    _check_parameter(J.m)
    return _jacobi_am(u, J.landen)
end
