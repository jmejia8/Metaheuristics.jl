"""
    PolynomialMutation(;η, p, bounds)

Polynomial mutation.
"""
struct PolynomialMutation
    η::Float64
    p::Float64
    bounds
    rng
end

function PolynomialMutation(;
        η = 15.0,
        bounds = zeros(0,0),
        p=1/getdim(bounds),
        rng = default_rng_mh()
    )
    PolynomialMutation(η, isfinite(p) ? p : 1e-2, bounds, rng)
end

function mutation!(Q, parameters::PolynomialMutation)
    bounds = parameters.bounds
    η = parameters.η
    p = parameters.p
    rng = parameters.rng
    for i in 1:size(Q,1)
        polynomial_mutation!(view(Q, i,:), bounds, η, p, rng)
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
        η=20,
        prob = 1 / length(vector),
        rng = default_rng_mh()
    )

    do_mutation = rand(rng, length(vector)) .< prob

    xu = view(bounds.ub, do_mutation)
    xl = view(bounds.lb, do_mutation)
    x = view(vector, do_mutation)

    δ1 = (x - xl) ./ (xu - xl)
    δ2 = (xu - x) ./ (xu - xl)

    D = length(xu)
    R = rand(rng, D)
    mask = R .< 0.5
    s = η+1
    mut_pow = 1 / (η + 1)
    δq = [ mask[i] ?
            ^(2R[i] + (1. - 2R[i]) * ^(1 - δ1[i], s), mut_pow) - 1 :
            1 - (2 * (1 - R[i]) + 2 * (R[i] - 0.5) * ^(1 - δ2[i], s))^mut_pow
            for i in 1:D
        ]

    vector[do_mutation] = _to_int_if_necessary(eltype(vector), x + δq .* ( xu - xl))
    vector
end

function polynomial_mutation!(vector,
        bounds::AbstractMatrix,
        η=20,
        prob = 1 / length(vector),
        rng = default_rng_mh()
    )

    b = BoxConstrainedSpace(lb = bounds[1,:], ub = bounds[2,:])
    polynomial_mutation!(vector, b, η, prob, rng)
end

