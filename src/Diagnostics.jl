# Self-checks for constructed orbits: every trajectory is compared with the geodesic
# equations by finite differences (independently of the member's own velocity and
# residual functions), the polar azimuth over one polar period with a quadrature of its own,
# and evaluation is timed. Used by `kerr_geo_diagnose`.

"""Mino-time rates (dt/dλ, dφ/dλ) of a Kerr geodesic at (r, z = cos θ), with sin²θ = 1 − z²."""
function _geodesic_rates(a, energy, lz, r, z)
    Δ = kerr_delta(a, r)
    P = kerr_radial_momentum(a, energy, lz, r)
    s2 = 1 - z^2
    dt = (r^2 + a^2) * P / Δ + a * lz - a^2 * energy * s2
    dphi = a * P / Δ - a * energy + (iszero(lz) ? 0.0 : lz / s2)
    return dt, dphi, s2
end

_first_field(nt, names) = begin
    for n in names
        hasproperty(nt, n) || continue
        f = getproperty(nt, n)
        f isa Function && return f
    end
    nothing
end

# Mino-time window on which to sample a member (open ends are pulled in by `margin`).
function _diagnostic_domain(m)
    dom = nothing
    if hasproperty(m, :Domain) && m.Domain isa NamedTuple && haskey(m.Domain, :mino)
        dom = m.Domain.mino
    elseif hasproperty(m, :Status) && haskey(m.Status, :domain) && m.Status.domain isa NamedTuple &&
            haskey(m.Status.domain, :mino) && m.Status.domain.mino isa Tuple
        dom = m.Status.domain.mino
    elseif hasproperty(m, :Status) && haskey(m.Status, :duration)
        dom = (0.0, m.Status.duration.horizon_lambda)
    end
    return dom === nothing ? (0.0, Inf) : float.(dom)
end

function _diagnostic_window(m; span=12.0, margin=0.03)
    lo, hi = _diagnostic_domain(m)
    lo = isfinite(lo) ? lo : (isfinite(hi) ? hi - span : -span / 2)
    hi = isfinite(hi) ? hi : lo + span
    w = hi - lo
    return (lo + margin * w, hi - margin * w)
end

# Far Mino times on the infinite sides of the domain: long-time and asymptotic regimes
# (repeated-root approach, large radii) where closed forms tend to lose precision.
function _diagnostic_far(m; offsets=(30.0, 60.0, 120.0))
    lo, hi = _diagnostic_domain(m)
    wlo, whi = _diagnostic_window(m)
    far = Float64[]
    isfinite(lo) || append!(far, whi .- offsets)
    isfinite(hi) || append!(far, wlo .+ offsets)
    return far
end

# Samples on the way to an infinity end reached at finite λ∞: r = 1e6, 1e12, 1e18, 1e24,
# located by bisection on r(λ), with a stencil scaled to λ∞ − λ. Returns (samples, radii
# that λ cannot resolve any more).
function _diagnostic_tail(m, rf)
    lo, hi = _diagnostic_domain(m)
    roles = hasproperty(m, :Domain) && haskey(m.Domain, :endpoint_roles) ?
        m.Domain.endpoint_roles : ()
    samples = Tuple{Float64,Float64}[]
    unresolved = 0
    for (end_λ, role, inner) in ((lo, isempty(roles) ? :none : roles[1], 1),
            (hi, isempty(roles) ? :none : roles[2], -1))
        (isfinite(end_λ) && occursin("infinity", String(role))) || continue
        wlo, whi = _diagnostic_window(m)
        start = inner > 0 ? whi : wlo                    # an interior point on that side
        for target in (1e6, 1e12, 1e18, 1e24)
            a_, b_ = start, end_λ
            r0 = try rf(a_) catch; continue end
            r0 < target || continue
            for _ in 1:400
                mid = (a_ + b_) / 2
                (mid == a_ || mid == b_) && break
                rm = try rf(mid) catch; Inf end
                (isfinite(rm) && rm < target) ? (a_ = mid) : (b_ = mid)
            end
            h = 1e-3 * abs(end_λ - a_)
            # λ no longer resolves r there, or the finite-difference nodes λ ± kh fall
            # within a few thousand ulps of each other
            if (try rf(a_) < target / 2 catch; true end) || h < 1e3 * eps(a_)
                unresolved += 1
            else
                push!(samples, (a_, h))
            end
        end
    end
    return samples, unresolved
end

_fd4(f, λ, h) = (-f(λ + 2h) + 8f(λ + h) - 8f(λ - h) + f(λ - 2h)) / (12h)
# Five-point derivative with a Richardson error estimate: (value, estimated error).
function _fd(f, λ, h)
    d1 = _fd4(f, λ, h); d2 = _fd4(f, λ, h / 2)
    d = d2 + (d2 - d1) / 15
    # the nodes λ ± kh are rounded to ulp(λ): a relative error of ~ulp(λ)/h in the step,
    # which dominates far out on a finite Mino-time domain (λ close to λ∞)
    return (d, abs(d2 - d1) + 8eps(λ) / h * abs(d))
end

"""
    _polar_period_azimuth(a, E, Lz, Q, sector)

Δφ of the polar motion over one period of z² (one pass between the pendular turning points,
two for vortical motion), by its own quadrature: with ζ measured from the turning point
nearest the axis, sin²θ = ε + u sin²ζ keeps the digits of ε = 1 − z²_turn and the square
root of Θ cancels against dz, so the integrand is smooth up to the Lz/ε spike, which
the adaptive rule resolves. For a = 0 the pendular value is π sgn(Lz) exactly.
"""
function _polar_period_azimuth(a, energy, lz, q, sector)
    T = _float_type(a, energy, lz, q)
    quarter, rtol = T(π) / 2, _tol(T, 1e-14)
    c = -a^2 * _e2m1(energy)                          # Θ(u) = c u² − (Q + Lz² + c) u + Q
    b = q + lz^2 + c
    if sector === :vortical
        # u₋ < u₊ < 1, c < 0; 1 − u₊ = −Lz²/(c (1 − u₋)) from Θ(1) = −Lz²
        # u₊ from the stable pair of Θ, 1 − u₊ from that of Θ(1 − w) = c w² − (2c − b) w − Lz²
        # (same discriminant) and the width u₊ − u₋ = s/|c| directly: each keeps its digits whether
        # the band lies next to the equator or next to the axis
        s = sqrt(b^2 - 4c * q)
        big_root = (b + copysign(s, b)) / (2c)
        up = max(big_root, q / (c * big_root))
        t = 2c - b
        ε = -lz^2 / (c * ((t + copysign(s, t)) / (2c)))    # 1 − u₊ = −Lz²/(c (1 − u₋))
        dw = s / (-c)
        β = -c
        f(ζ) = lz / (sqrt(up - dw * sin(ζ)^2) * (ε + dw * sin(ζ)^2) * sqrt(β))
        return 2 * quadgk(f, zero(T), quarter; rtol=rtol)[1]
    end
    iszero(a) && return T(π) * sign(lz)
    # pendular: u₁ = z²_max is the root of Θ in (0, 1]; G(u) = Θ(u)/(u₁ − u) = Q/u₁ − c u.
    # ε = 1 − u₁ is the root in (0, 1] of Θ(1 − w) = c w² − (2c − b) w − Lz² (same discriminant,
    # roots with product −Lz²/c), so it keeps its digits when u₂ lies within rounding of 1
    if iszero(c)
        u1 = q / b; ε = lz^2 / b
    else
        s = sqrt(b^2 - 4c * q)
        r_big = (b + copysign(s, b)) / (2c)
        r_small = q / (c * r_big)
        u1 = 0 < r_small <= 1 ? r_small : r_big
        t = 2c - b
        w_big = (t + copysign(s, t)) / (2c)
        w_small = -lz^2 / (c * w_big)
        ε = 0 < w_small <= 1 ? w_small : w_big
    end
    g(ζ) = lz / ((ε + u1 * sin(ζ)^2) * sqrt(q / u1 - c * u1 * cos(ζ)^2))
    return 2 * quadgk(g, zero(T), quarter; rtol=rtol)[1]
end

"""
    kerr_geo_diagnose(member; a, E, Lz, Q, samples=41, tol=1e-6, slow_eval=1e-3, far_samples=true)
    kerr_geo_diagnose(family::KerrGeodesicFamily; kwargs...)

Check an orbit against the geodesic equations using only its own callables, and time its
evaluation. On `samples` equally spaced Mino times across the member's domain (a stretch of
length 12 where the domain is infinite) and, when `far_samples`, three more on each infinite
side (30, 60 and 120 from the opposite edge of that stretch) and at r = 1e6, 1e12, 1e18,
1e24 on the way to an infinity end reached at finite λ, it compares, by finite differences
of the trajectory,

- `(dr/dλ)² = R(r)` and `(dz/dλ)² = Θ(z)` (relative to the size of the terms),
- `dt/dλ`, `dφ/dλ` with the Carter equations (or `dv/dλ`, `dψ/dλ` for members given in
  ingoing coordinates).

Each check reports the largest raw relative discrepancy over the samples it could decide
(`radial`, `polar`, `time`, `azimuth`) and, in `allowance`, the largest margin granted
there: the finite-difference error (Richardson estimate), the `eps·|t|` resolution of the
values, and the precision of the reference rates next to a horizon (r − r₊) or next to the
axis (1 − z²). A sample passes when discrepancy ≤ `tol` + allowance. Samples whose
allowance alone exceeds `tol` cannot decide the check (the stencil does not resolve the
rate, e.g. the Lz/sin²θ spike of a near-axis orbit) and are counted in `skipped`, as are
tail radii that λ cannot resolve (`skipped.tail`). `period_azimuth` compares the polar Δφ
over one polar period with an independent quadrature (NaN when the polar motion does not
oscillate or Lz = 0). It also checks that `t` (or `v`) increases, that `v − t − r*` and
`ψ − φ − φ_H` stay constant when both charts exist, and counts non-finite values.

Returns a NamedTuple: `member` (the case ID), `ok` (no issue other than slow evaluation),
`slow` (median single evaluation above `slow_eval` seconds), `issues`, the discrepancies
`radial`, `polar`, `time`, `azimuth` with their `allowance`, `tolerance`, `skipped`,
`ingoing_consistency` (largest relative drift of `v − t − r*` and `ψ − φ − φ_H`; 0 with a
single chart), `period_azimuth`, the number of checked `samples`, `eval_median` and
`eval_max`. A family returns one record per member; members that could not be built are
listed in `family.Status.member_errors`.
"""
function kerr_geo_diagnose(m; a, E, Lz, Q, samples::Int=41, tol=1e-6,
        slow_eval=1e-3, far_samples::Bool=true)
    tr = hasproperty(m, :Trajectory) ? m.Trajectory : nothing
    name = hasproperty(m, :CaseId) ? m.CaseId :
           hasproperty(m, :Formula) ? m.Formula : nameof(typeof(m))
    issues = String[]
    empty = (radial=NaN, polar=NaN, time=NaN, azimuth=NaN)
    tr === nothing && return (member=name, ok=false, slow=false, issues=["no trajectory"],
        empty..., ingoing_consistency=NaN, period_azimuth=NaN, tolerance=tol,
        allowance=empty, skipped=(radial=0, polar=0, time=0, azimuth=0, tail=0),
        samples=0, eval_median=NaN, eval_max=NaN)
    rf = _first_field(tr, (:r,))
    zf = _first_field(tr, (:z, :q))
    θf = _first_field(tr, (:theta, :θ))
    zf === nothing && θf !== nothing && (zf = λ -> cos(θf(λ)))
    tf = _first_field(tr, (:t,)); φf = _first_field(tr, (:phi, :ϕ))
    vf = _first_field(tr, (:v,)); ψf = _first_field(tr, (:psi,))
    ingoing = false
    if (tf === nothing || φf === nothing) && vf !== nothing && ψf !== nothing
        tf, φf, ingoing = vf, ψf, true
    end
    lo, hi = _diagnostic_window(m)
    h0 = 1e-4 * max(1.0, hi - lo)
    pts = [(λ, h0) for λ in range(lo, hi; length=samples)]
    far_samples && append!(pts, [(λ, h0) for λ in _diagnostic_far(m)])
    tail, tail_unresolved = far_samples ? _diagnostic_tail(m, rf) : (Tuple{Float64,Float64}[], 0)
    append!(pts, tail)
    sort!(pts; by=first)
    inside(λ, h) = λ - 2h >= lo - 0.03(hi - lo) && λ + 2h <= hi + 0.03(hi - lo)
    tail_λ = Set(first.(tail)); far_λ = Set(_diagnostic_far(m))
    both_charts = !ingoing && vf !== nothing && ψf !== nothing
    charts = Tuple{NTuple{2,Float64},Float64}[]   # (v − t − r*, ψ − φ − φ_H) and its scale
    rp = _rplus(a)
    errs = Dict(k => 0.0 for k in (:radial, :polar, :time, :azimuth))
    allow = Dict(k => 0.0 for k in (:radial, :polar, :time, :azimuth))
    skipped = Dict(k => 0 for k in (:radial, :polar, :time, :azimuth))
    failed = Dict(k => 0.0 for k in (:radial, :polar, :time, :azimuth))
    # one check: raw discrepancy `mis` against `allowance` (absolute), relative to `scale`.
    # A sample whose allowance is above the tolerance or not finite (a nonfinite reference
    # rate) cannot decide the check and is counted as skipped; a nonfinite discrepancy
    # against a finite allowance is a failure.
    function record!(k, mis, allowance, scale)
        raw = mis / scale; al = allowance / scale
        if !(al <= tol)
            skipped[k] += 1
            return
        end
        isfinite(raw) || (raw = Inf)
        errs[k] = max(errs[k], raw); allow[k] = max(allow[k], al)
        raw > tol + al && (failed[k] = max(failed[k], raw))
    end
    times = Float64[]
    lastt = -Inf; monotone = true; nonfinite = 0; nonfinite_t = 0; checked = 0
    for (λ, h) in pts
        (λ in tail_λ || λ in far_λ || inside(λ, h)) || continue
        local r, z, dr, dz, edr, edz
        try
            t0 = time_ns(); r = rf(λ); push!(times, (time_ns() - t0) * 1e-9)
            z = zf === nothing ? 0.0 : zf(λ)
            (dr, edr) = _fd(rf, λ, h)
            (dz, edz) = zf === nothing ? (0.0, 0.0) : _fd(zf, λ, h)
        catch err
            push!(issues, "λ=$(round(λ; sigdigits=4)): " * first(split(sprint(showerror, err), '\n')))
            continue
        end
        if !all(isfinite, (r, z, dr, dz))
            nonfinite += 1; continue
        end
        checked += 1
        R = kerr_radial_potential(a, E, Lz, Q, r)
        Θ = kerr_polar_z_potential(a, E, Lz, Q, z)
        record!(:radial, abs(dr^2 - R), 100abs(dr) * edr, max(1.0, r^4))
        record!(:polar, abs(dz^2 - Θ), 100abs(dz) * edz, max(1.0, abs(Q) + Lz^2 + a^2))
        if tf !== nothing && φf !== nothing && r > rp * (1 + 1e-6)
            try
                t0 = time_ns(); tv = tf(λ); push!(times, (time_ns() - t0) * 1e-9)
                (dt, edt) = _fd(tf, λ, h); (dφ, edφ) = _fd(φf, λ, h)
                Tt, Tφ, s2 = _geodesic_rates(a, E, Lz, r, z)
                # a value of size |t| is known to eps|t| at best: finite differences of it
                # cannot resolve rates below ~eps|t|/h (large t near apoapsis/infinity)
                edt += 8eps() * abs(tv) / h; edφ += 8eps() * abs(φf(λ)) / h
                if ingoing                     # v = t + r*, ψ = φ + horizon azimuth
                    Δ = kerr_delta(a, r)
                    Tt += (r^2 + a^2) / Δ * dr; Tφ += a / Δ * dr
                    edt += abs((r^2 + a^2) / Δ) * edr; edφ += abs(a / Δ) * edr
                end
                all(isfinite, (tv, dt, dφ)) || (nonfinite_t += 1)
                # the reference rates are only as good as r − r₊ next to a horizon, and the
                # polar φ rate Lz/(1 − z²) only as good as 1 − z² from a rounded z
                refrel = 8eps() * max(1.0, r) / max(r - rp, eps())
                polar_ref = iszero(Lz) ? 0.0 : 4eps() / max(s2, floatmin()) * abs(Lz / s2)
                if all(isfinite, (dt, dφ))
                    record!(:time, abs(dt - Tt), 50edt + refrel * abs(Tt), max(1.0, abs(Tt)))
                    # on the axis (Lz = 0, sin θ = 0) φ is a gauge, not a coordinate of the orbit
                    if iszero(Lz) && s2 < 1e-10
                        skipped[:azimuth] += 1
                    else
                        record!(:azimuth, abs(dφ - Tφ), 50edφ + refrel * abs(Tφ) + polar_ref,
                            max(1.0, abs(Tφ)))
                    end
                end
                if both_charts       # v − t − r*, ψ − φ − φ_H are constants of the motion
                    # (v, ψ may be defined on part of the domain only, e.g. a Class N member's
                    # future half; on the axis φ and ψ are gauges and only v is checked)
                    cv = try
                        (vf(λ) - tv - kerr_rstar(a, r),
                         iszero(Lz) && s2 < 1e-10 ? 0.0 : ψf(λ) - φf(λ) - _horizon_azimuth(a, r))
                    catch err
                        err isa DomainError || rethrow()
                        (NaN, NaN)
                    end
                    all(isfinite, cv) && push!(charts, (cv, max(1.0, abs(tv))))
                end
                isfinite(tv) && (tv < lastt - 1e-9 * max(1, abs(tv)) && (monotone = false); lastt = tv)
            catch err
                push!(issues, "t/φ at λ=$(round(λ; sigdigits=4)): " *
                    first(split(sprint(showerror, err), '\n')))
            end
        end
    end
    # the polar azimuth over one polar period against an independent quadrature
    period_azimuth = NaN
    meta = hasproperty(m, :Status) && haskey(m.Status, :polar) ? m.Status.polar : nothing
    if meta isa NamedTuple && haskey(meta, :mean_rates) && haskey(meta, :period_u) &&
            !iszero(Lz) && get(meta, :sector, :none) !== :axis_crossing
        sector = meta.sector === :vortical ? :vortical : :pendular
        got = meta.mean_rates.phi * meta.period_u / meta.omega
        ref = _polar_period_azimuth(a, E, Lz, Q, sector)
        period_azimuth = abs(got - ref) / max(1.0, abs(ref))
        period_azimuth > tol && push!(issues,
            "polar Δφ per period off by $(round(period_azimuth; sigdigits=3))")
    end
    nonfinite > 0 && push!(issues, "$nonfinite non-finite samples")
    nonfinite_t > 0 && push!(issues, "$nonfinite_t non-finite t/φ samples")
    monotone || push!(issues, (ingoing ? "v" : "t") * " decreases")
    # compared with the sample where the constants are known best (smallest |t|): far out
    # t and r* are large and v − t − r* keeps only eps·|t| of them
    consist = 0.0
    if !isempty(charts)
        cref, sref = charts[argmin(last.(charts))]
        consist = maximum(maximum(abs.(cv .- cref)) / max(sc, sref) for (cv, sc) in charts)
    end
    consist > 100tol && push!(issues, "ingoing/BL mismatch $(round(consist; sigdigits=3))")
    for k in (:radial, :polar, :time, :azimuth)
        failed[k] > 0 && push!(issues, "$k equation error $(round(failed[k]; sigdigits=3))" *
            " (allowance $(round(allow[k]; sigdigits=3)))")
    end
    emax = isempty(times) ? NaN : maximum(times)
    emed = isempty(times) ? NaN : sort(times)[cld(length(times), 2)]
    slow = !isempty(times) && emed > slow_eval
    slow && push!(issues, "slow evaluation: median $(round(emed * 1e6; sigdigits=3)) μs")
    nt(d) = (radial=d[:radial], polar=d[:polar], time=d[:time], azimuth=d[:azimuth])
    return (member=name, ok=isempty(filter(s -> !startswith(s, "slow"), issues)),
            slow=slow, issues=issues, nt(errs)..., ingoing_consistency=consist,
            period_azimuth=period_azimuth, tolerance=tol, allowance=nt(allow),
            skipped=merge(nt(skipped), (tail=tail_unresolved,)), samples=checked,
            eval_median=emed, eval_max=emax)
end

function kerr_geo_diagnose(f::KerrGeodesicFamily; kwargs...)
    a = f.Parameters.a
    E, Lz, Q = f.ConstantsOfMotion.E, f.ConstantsOfMotion.Lz, f.ConstantsOfMotion.Q
    return [kerr_geo_diagnose(m; a=a, E=E, Lz=Lz, Q=Q, kwargs...) for m in kerr_geo_members(f)]
end
