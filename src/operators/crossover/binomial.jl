struct BinomialCrossover{R}
    p::Float64
    n_offsprings::Int
    rng::R
end

BinomialCrossover(p, n_offsprings, rng) = BinomialCrossover{typeof(rng)}(
    p,
    n_offsprings,
    rng,
)

function BinomialCrossover(; p = 0.5, n_offsprings = 2, rng = default_rng_mh())
    BinomialCrossover(p, n_offsprings, rng)
end

function crossover(population, c::BinomialCrossover)
    Q = positions(population)
    n = size(Q, 1) ÷ 2
    R = rand(c.rng, n, size(Q, 2))
    p = c.p

    if c.n_offsprings > 2
        # odd rows, taking the genes of the even rows where R .<= p
        Q1 = Q[1:2:end - 1, :]
        @inbounds for j in axes(R, 2), i in 1:n
            Q1[i, j] = ifelse(R[i, j] <= p, Q[2i, j], Q1[i, j])
        end
        return Q1
    end

    # first half: odd rows (unchanged); second half: even rows, taking the genes of
    # the odd rows where R .<= p
    O = similar(Q, 2n, size(Q, 2))
    @inbounds for j in axes(R, 2), i in 1:n
        O[i, j] = Q[2i - 1, j]
        O[n + i, j] = ifelse(R[i, j] <= p, Q[2i - 1, j], Q[2i, j])
    end
    return O
end
