"""
    BinomialCrossover(;p = 0.5, n_offsprings = 2)

Binomial crossover between consecutive pairs of parents (rows `2i-1` and `2i`). The first
child of each pair takes the gene of the first parent with probability `p`, otherwise the
gene of the second parent. With `n_offsprings = 1` only that child is returned (one per
pair, as in BRKGA). With `n_offsprings = 2` its complement is returned too: the first `n`
rows are the first children and the last `n` rows the second children.
"""
struct BinomialCrossover{R}
    p::Float64
    n_offsprings::Int
    rng::R

    function BinomialCrossover{R}(p, n_offsprings, rng::R) where {R}
        n_offsprings in (1, 2) ||
            throw(ArgumentError("n_offsprings must be 1 or 2, got $n_offsprings"))
        new{R}(p, n_offsprings, rng)
    end
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

    # first children: gene of the first parent (odd row) where R .<= p, else of the
    # second parent (even row); second children (n_offsprings = 2): the complement
    O = similar(Q, c.n_offsprings * n, size(Q, 2))
    if c.n_offsprings == 1
        @inbounds for j in axes(R, 2), i in 1:n
            O[i, j] = ifelse(R[i, j] <= p, Q[2i - 1, j], Q[2i, j])
        end
    else
        @inbounds for j in axes(R, 2), i in 1:n
            first = R[i, j] <= p
            a, b = Q[2i - 1, j], Q[2i, j]
            O[i, j] = ifelse(first, a, b)
            O[n + i, j] = ifelse(first, b, a)
        end
    end
    return O
end
