"""
    SBX(;η, p, bounds)

Simulated Binomial Crossover.
"""
struct SBX{B, R}
    η::Float64
    p::Float64
    bounds::B
    rng::R
end

SBX(η, p, bounds, rng) = SBX{typeof(bounds), typeof(rng)}(η, p, bounds, rng)

SBX(; η = 15, p = 0.9, bounds = zeros(0, 0), rng = default_rng_mh()) = SBX(
    η,
    p,
    _mat_to_bounds(bounds),
    rng,
)

function crossover(population, parameters::SBX)
    isempty(population) && return zeros(0, 0)
    Q = positions(population)
    bounds = parameters.bounds
    rng = parameters.rng
    for i in 1:2:length(population) - 1
        # rows i and i+1 of Q hold the parents and are replaced by the offspring
        c1 = view(Q, i, :)
        c2 = view(Q, i + 1, :)
        SBX_crossover!(c1, c2, bounds, parameters.η, parameters.p, rng)
    end
    Q
end

function _sbx_βq(β, η, R)
    α = 2 - β^(-η - 1)
    s = 1 / (η + 1)
    return R <= 1 / α ? (R * α)^s : (1 / (2 - R * α))^s
end

"""
    SBX_crossover(vector1, vector2, bounds, η=15, p_variable = 0.9)

Simulated binomial crossover for given two `Vectors{Real}`.
"""
function SBX_crossover(
    vector1,
    vector2,
    bounds::BoxConstrainedSpace,
    η = 15,
    p_variable = 0.9,
    rng = default_rng_mh(),
)
    SBX_crossover!(copy(vector1), copy(vector2), bounds, η, p_variable, rng)
end

"""
    SBX_crossover!(c1, c2, bounds, η=15, p_variable = 0.9)

In-place version of [`SBX_crossover`](@ref): `c1` and `c2` contain the parents on input and
the offspring on output.
"""
function SBX_crossover!(
    c1,
    c2,
    bounds::BoxConstrainedSpace,
    η = 15,
    p_variable = 0.9,
    rng = default_rng_mh(),
)

    xu = bounds.ub
    xl = bounds.lb
    D = getdim(bounds)
    if length(c1) != D || length(c2) != D
        throw(DimensionMismatch("parents and bounds must have the same dimension $D"))
    end

    # Random numbers are drawn in bulk (same order as previous versions), so seeded runs
    # are reproducible across versions.
    u = rand(rng, D)
    R = rand(rng, D)
    swap = rand(rng, Bool, D)

    @inbounds for j in 1:D
        a = c1[j]
        b = c2[j]
        (u[j] > p_variable || abs(b - a) <= eps()) && continue

        y1 = min(a, b)
        y2 = max(a, b)
        Δ = max(eps(), y2 - y1)

        # bounded SBX (Deb & Agrawal, 1995): child 1 depends on the distance to the lower
        # bound and child 2 on the distance to the upper bound
        βq = _sbx_βq(1 + (2 * (y1 - xl[j]) / Δ), η, R[j])
        o1 = 0.5 * (y1 + y2 - βq * Δ)

        βq = _sbx_βq(1 + (2 * (xu[j] - y2) / Δ), η, R[j])
        o2 = 0.5 * (y1 + y2 + βq * Δ)

        if swap[j]
            o1, o2 = o2, o1
        end

        c1[j] = _to_int_if_necessary(eltype(c1), o1)
        c2[j] = _to_int_if_necessary(eltype(c2), o2)
    end

    reset_to_violated_bounds!(c1, bounds)
    reset_to_violated_bounds!(c2, bounds)

    return c1, c2
end

function SBX_crossover(
    vector1,
    vector2,
    bounds::AbstractMatrix,
    η = 15,
    p_variable = 0.9,
    rng = default_rng_mh(),
)
    b = BoxConstrainedSpace(lb = bounds[1, :], ub = bounds[2, :])
    SBX_crossover(vector1, vector2, b, η, p_variable, rng)
end
