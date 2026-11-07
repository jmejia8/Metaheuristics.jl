"""
    UniformCrossover(;p = 0.5)

Uniform crossover a.k.a. Binomial crossover. Suitable for binary representation.
"""
struct UniformCrossover{R}
    p::Float64
    rng::R
end

UniformCrossover(p, rng) = UniformCrossover{typeof(rng)}(p, rng)

UniformCrossover(; p = 0.5, rng = default_rng_mh()) = UniformCrossover(p, rng)

function crossover(population, parameters::UniformCrossover)
    n = length(population) ÷ 2
    rng = parameters.rng
    offspring_A = positions(population[1:n])
    offspring_B = positions(population[n + 1:2n])
    mask = rand(rng, size(offspring_A)...) .<= parameters.p
    # swapping genes
    tmp = copy(offspring_A[mask])
    offspring_A[mask] = offspring_B[mask]
    offspring_B[mask] = tmp
    [offspring_A; offspring_B]
end

"""
    DE_crossover(x, u, CR)

Binomial crossover between x and u for Differential Evolution with probability CR, i.e.,
`v[j] = u[j]` if `rand() < CR`, otherwise `v[j] = x[j]`. Return `v`.
"""
function DE_crossover(x, u, CR)
    D = length(x)
    # binomial crossover
    v = zeros(D)
    j_rand = rand(1:D)

    # binomial crossover
    for j in 1:D
        if rand() < CR || j == j_rand
            v[j] = u[j]
        else
            v[j] = x[j]
        end
    end

    return v
end
