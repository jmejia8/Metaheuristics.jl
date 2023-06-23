mutable struct BinomialCrossover
    p::Float64
    n_offsprings::Int
    rng
end

function BinomialCrossover(;p = 0.5, n_offsprings = 2, rng=default_rng_mh())
    BinomialCrossover(p, n_offsprings, rng)
end

function crossover(population, c::BinomialCrossover)
    Q = positions(population)

    Q1 = Q[1:2:end-1, :]
    Q2 = Q[2:2:end,   :]

    mask = rand(c.rng, size(Q,1)÷2, size(Q, 2)) .<= c.p
    Q1[mask] = Q2[mask]

    if c.n_offsprings > 2
        return Q1
    end

    Q1 = Q[1:2:end-1, :]
    Q2[mask] = Q1[mask]
    return vcat(Q1, Q2)
end

