"""
    BitFlipMutation(;p = 1e-2)

Flip each bit with probability `p`.
"""
struct BitFlipMutation
    p::Float64
    BitFlipMutation(;p = 1e-2) = new(p)
end

function mutation!(Q, parameters::BitFlipMutation)
    mask = rand(size(Q)...) .< parameters.p
    Q[mask] = .!Q[mask]
    Q
end

"""
    SlightMutation

> Fogel, D. B. (1988). An evolutionary approach to the traveling salesman problem.
> Biological Cybernetics, 60(2), 139-144.
"""
struct SlightMutation end

function mutation!(Q, parameters::SlightMutation)
    N, D = size(Q)
    k = rand(1:D, N)
    s = rand(1:D, N)
    for i in 1:N
        s[i] == k[i] && continue
        if s[i] < k[i]
            c = vcat(1:s[i]-1, k[i], s[i]:k[i]-1, k[i]+1:D)
        else s[i] > k[i]
            c = vcat(1:k[i]-1, k[i]+1:s[i]-1, k[i], s[i]:D)
        end
        Q[i,:] = Q[i, c]
    end
    Q
end


######################

"""
    polynomial_mutation!(vector, bounds, η=20, prob = 1 / length(vector))

Polynomial Mutation applied to a vector of real numbers.
"""
function polynomial_mutation!(vector, bounds::Bounds, η=20, prob = 1 / length(vector))
    do_mutation = rand(length(vector)) .< prob

    xu = view(bounds.ub, do_mutation)
    xl = view(bounds.lb, do_mutation)
    x = view(vector, do_mutation)

    δ1 = (x - xl) ./ (xu - xl)
    δ2 = (xu - x) ./ (xu - xl)

    D = length(xu)
    R = rand(D)
    mask = R .< 0.5
    s = η+1.0
    mut_pow = 1.0 / (η + 1.0)
    δq = [ mask[i] ?
            ^(2.0R[i] + (1. - 2.0R[i]) * ^(1.0 - δ1[i], s), mut_pow) - 1.0 :
            1.0 - (2.0 * (1.0 - R[i]) + 2.0 * (R[i] - 0.5) * ^(1.0 - δ2[i], s))^mut_pow
            for i in 1:D
        ]

    vector[do_mutation] = _to_int_if_necessary(eltype(vector), x + δq .* ( xu - xl))
    vector
end

function polynomial_mutation!(vector, bounds::AbstractMatrix, η=20, prob = 1 / length(vector))
    b = Bounds(lb = bounds[1,:], ub = bounds[2,:])
    polynomial_mutation!(vector, b, η, prob)
end

"""
    PolynomialMutation(;η, p, bounds)

Polynomial mutation.
"""
struct PolynomialMutation
    η::Float64
    p::Float64
    bounds
    PolynomialMutation(;η = 15.0, bounds = zeros(0,0), p=1/getdim(bounds)) = begin
        new(η, isfinite(p) ? p : 1e-2, bounds)
    end
end

function mutation!(Q, parameters::PolynomialMutation)
    for i in 1:size(Q,1)
        polynomial_mutation!(view(Q, i,:), parameters.bounds, parameters.η, parameters.p)
    end
    Q
end

"""
    DE_mutation(population, F = 1.0, strategy = :rand1)

Generate a `Vector` computed from population used in Differential Evolution.
Parameters: `F` is the stepsize, `strategy` can be one the following `:best1`, `:rand2`,
`:randToBest1` or `:best2`.
"""
function DE_mutation(population,
                    F::Float64 = 1.0,
                    strategy::Symbol=:rand1,
                    best_ind=0
    )

    N = length(population)

    i = rand(1:N)
    # select participats
    r1 = rand(1:N)
    while r1 == i
        r1 = rand(1:N)
    end

    r2 = rand(1:N)
    while r2 == i || r1 == r2
        r2 = rand(1:N)
    end


    x = get_position(population[i] )
    a = get_position(population[r1])
    b = get_position(population[r2])

    # strategy is selected here
    if strategy == :rand1
        # DE/rand/1
        return  x + F * (a - b)
    end

    r3 = rand(1:N)
    while r3 == i || r3 == r1 || r3 == r2
        r3 = rand(1:N)
    end

    c = get_position(population[r3])

    if strategy == :rand2
        # DE/rand/2

        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3
            r4 = rand(1:N)
        end

        r5 = rand(1:N)
        while r5 == i || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4
            r5 = rand(1:N)
        end

        d = get_position(population[r4])
        ee = get_position(population[r5])

        return ee + F * (a - b + c - d)
    end

    best_ind = best_ind == 0 ? argbest(population) : best_ind
    x_best = get_position(population[best_ind])
    if strategy == :best1
        # DE/best/1
        u = x_best + F * (b - c)
    elseif strategy == :randToBest1
        # DE/rand-to-best/1
        u = x + F * (x_best - x + a - b)
    elseif strategy == :best2
        # DE/best/2
        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3 || r4 == best_ind
            r4 = rand(1:N)
        end
        d = population[r4].x
        u = x_best + F * (a - b + c - d)
    else
        error("Unknown strategy " * string(strategy))
    end

    return u
end

