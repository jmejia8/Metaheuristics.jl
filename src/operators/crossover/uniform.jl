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
    # rows 1:n are the first parents, rows n+1:2n the second ones
    O = positions(view(population, 1:2n))
    R = rand(parameters.rng, n, size(O, 2))
    p = parameters.p
    # swapping genes (branchless: R is random, so a branch is mispredicted often)
    @inbounds for j in axes(R, 2), i in 1:n
        swap = R[i, j] <= p
        a, b = O[i, j], O[n + i, j]
        O[i, j] = ifelse(swap, b, a)
        O[n + i, j] = ifelse(swap, a, b)
    end
    O
end

"""
    DE_crossover(x, u, CR; rng)

Binomial crossover between x and u for Differential Evolution with probability CR, i.e.,
`v[j] = u[j]` if `rand() < CR`, otherwise `v[j] = x[j]`. Return `v`.
"""
function DE_crossover(x, u, CR; rng = default_rng_mh())
    v = zeros(length(x))
    v .= u
    DE_crossover!(v, x, CR; rng)
end

"""
    DE_crossover!(u, x, CR; rng)

In-place version of [`DE_crossover`](@ref): `u[j]` is replaced by `x[j]` for the genes
that are not crossed. Return `u`.
"""
function DE_crossover!(u, x, CR; rng = default_rng_mh())
    D = length(x)
    j_rand = rand(rng, 1:D)

    # binomial crossover
    @inbounds for j in 1:D
        if rand(rng) < CR || j == j_rand
            continue
        end
        u[j] = x[j]
    end

    return u
end
