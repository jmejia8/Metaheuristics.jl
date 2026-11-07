"""
    TournamentSelection(;K=2, N=0)

Perform the K-tournament selection and return N elements.
"""
struct TournamentSelection
    K::Int # K-tournament
    N::Int # N solutions
    TournamentSelection(; K = 2, N = 0) = new(K, N)
end

function selection(
    population,
    parameters::TournamentSelection;
    fitness = rank_solutions(population),
    rng = default_rng_mh(),
)
    parameters.K < 1 && throw(ArgumentError("TournamentSelection requires K ≥ 1."))
    n = length(population)
    out = Vector{Int}(undef, parameters.N)
    @inbounds for j in 1:parameters.N
        best = rand(rng, 1:n)
        for _ in 2:parameters.K
            i = rand(rng, 1:n)
            # same tie and NaN handling as argmin (first minimum wins)
            Base.isgreater(fitness[best], fitness[i]) && (best = i)
        end
        out[j] = best
    end
    out
end

"""
    RouletteWheelSelection(;N=0)

Perform Roulette Wheel Selection (minimization) and return N elements.

Solution `i` is chosen with probability proportional to `1 / fitness[i]`. If some fitness
value is ≤ 0, the values are shifted first so the smallest becomes 1, i.e., the weights are
`1 / (fitness[i] - minimum(fitness) + 1)`.
"""
struct RouletteWheelSelection
    N::Int # solutions required
    RouletteWheelSelection(; N = 0) = new(N)
end

function selection(
    population,
    parameters::RouletteWheelSelection;
    fitness = get_fitness(population),
    rng = default_rng_mh(),
)
    # positive weights: the smallest fitness maps to 1 when it is not positive
    fmin = minimum(fitness)
    shift = fmin > 0 ? eps() : fmin - 1
    p = 1 ./ (fitness .- shift)
    # cumulative sum, normalized
    @inbounds for i in 2:length(p)
        p[i] += p[i - 1]
    end
    p ./= maximum(p)
    # binary search requires a non-decreasing p without NaN
    if all(isfinite, p) && issorted(p)
        return [searchsortedlast(p, rand(rng)) + 1 for _ in 1:parameters.N]
    end
    [findfirst(>(rand(rng)), p) for _ in 1:parameters.N]
end

"""
    binary_tournament(population; rng)
Apply binary tournament to obtain a solution from from population.
"""
function binary_tournament(population; rng = default_rng_mh())
    a = rand(rng, 1:length(population))
    b = rand(rng, 1:length(population))

    return is_better(population[a], population[b]) ? population[a] : population[b]
end

function binary_tournament(
    population,
    fitness::Vector{T};
    rng = default_rng_mh(),
) where {T <: Real}
    a = rand(rng, 1:length(population))
    b = rand(rng, 1:length(population))

    return fitness[a] < fitness[b] ? population[a] : population[b]
end

mutable struct BiasedSelection
    num_elites::Int
    num_offsprings::Int
end

function selection(population, parameters::BiasedSelection)
    elites = 1:parameters.num_elites
    no_elites = parameters.num_elites + 1:length(population)
    # one (elite, non-elite) pair per offspring
    n = parameters.num_offsprings

    parent = ones(Int, 2n)
    parent[1:2:end] = rand(elites, n)
    parent[2:2:end] = rand(no_elites, n)
    parent
end
