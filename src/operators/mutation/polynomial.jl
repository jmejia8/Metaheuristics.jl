"""
    PolynomialMutation(;η, p, bounds)

Polynomial mutation.
"""
struct PolynomialMutation{B, R}
    η::Float64
    p::Float64
    bounds::B
    rng::R
end

PolynomialMutation(η, p, bounds, rng) = PolynomialMutation{typeof(bounds), typeof(rng)}(
    η,
    p,
    bounds,
    rng,
)

function PolynomialMutation(;
    η = 15.0,
    bounds = zeros(0, 0),
    p = 1 / getdim(_mat_to_bounds(bounds)),
    rng = default_rng_mh(),
)
    PolynomialMutation(η, isfinite(p) ? p : 1.0e-2, _mat_to_bounds(bounds), rng)
end

function mutation!(Q, parameters::PolynomialMutation)
    bounds = parameters.bounds
    η = parameters.η
    p = parameters.p
    rng = parameters.rng
    for i in 1:size(Q, 1)
        polynomial_mutation!(view(Q, i, :), bounds, η, p, rng)
    end
    Q
end

"""
    polynomial_mutation!(vector, bounds, η=20, prob = 1 / length(vector))

Polynomial Mutation applied to a vector of real numbers.
"""
function polynomial_mutation!(
    vector,
    bounds::BoxConstrainedSpace,
    η = 20,
    prob = 1 / length(vector),
    rng = default_rng_mh(),
)

    xu = bounds.ub
    xl = bounds.lb
    n = length(vector)
    if getdim(bounds) != n
        throw(DimensionMismatch("vector and bounds must have the same dimension $n"))
    end

    # Random numbers are drawn in bulk (same order as previous versions), so seeded runs
    # are reproducible across versions.
    u = rand(rng, n)
    R = rand(rng, count(<(prob), u))

    s = η + 1
    mut_pow = 1 / (η + 1)
    k = 0
    @inbounds for j in 1:n
        u[j] < prob || continue
        k += 1
        x = vector[j]
        Δ = xu[j] - xl[j]
        δ1 = (x - xl[j]) / Δ
        δ2 = (xu[j] - x) / Δ
        r = R[k]
        δq = if r < 0.5
            ^(2r + (1.0 - 2r) * ^(1 - δ1, s), mut_pow) - 1
        else
            1 - (2 * (1 - r) + 2 * (r - 0.5) * ^(1 - δ2, s))^mut_pow
        end
        vector[j] = _to_int_if_necessary(eltype(vector), x + δq * Δ)
    end
    vector
end

function polynomial_mutation!(
    vector,
    bounds::AbstractMatrix,
    η = 20,
    prob = 1 / length(vector),
    rng = default_rng_mh(),
)

    b = BoxConstrainedSpace(lb = bounds[1, :], ub = bounds[2, :])
    polynomial_mutation!(vector, b, η, prob, rng)
end
