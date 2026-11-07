"""
    OrderCrossover()
Order crossover for representation where order is important. Suitable for permutation
representation.
"""
struct OrderCrossover{R}
    rng::R
end

OrderCrossover(; rng = default_rng_mh()) = OrderCrossover(rng)

function crossover(population, parameters::OrderCrossover)
    O = positions(population)
    N, D = size(O)
    s = rand(parameters.rng, 1:D, N) # slash points
    for i in 1:2:N
        PA = O[i, :]     # parent A
        PB = O[i + 1, :] # parent B
        O[i, s[i] + 1:D] = setdiff(PB, PA[1:s[i]])
        O[i + 1, s[i] + 1:D] = setdiff(PA, PB[1:s[i]])
    end
    O
end
