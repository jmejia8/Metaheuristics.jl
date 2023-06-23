"""
    TournamentSelection(;K=2, N=0)

Perform the K-tournament selection and return N elements.
"""
struct TournamentSelection
    K::Int # K-tournament
    N::Int # N solutions
    TournamentSelection(;K=2, N=0) = new(K, N)
end

function selection(
        population,
        parameters::TournamentSelection;
        fitness = rank_solutions(population)
    )
    mask = 1:length(population)
    _winner(subset) = subset[argmin(fitness[subset])]
    [_winner(rand(mask, parameters.K)) for _ in 1:parameters.N]
end


"""
    RouletteWheelSelection(;N=0)

Perform Roulette Wheel Selection and return N elements.
"""
struct RouletteWheelSelection
    N::Int # solutions required
    RouletteWheelSelection(;N=0) = new(N)
end

function selection(
        population,
        parameters::RouletteWheelSelection;
        fitness=get_fitness(population)
    )
    fitness = fitness .- (min(minimum(fitness),zero(eltype(fitness))) + eps())
    fitness = cumsum(1 ./ fitness)
    p = fitness ./ maximum(fitness)
    [findfirst(>(rand()), p) for i in 1:parameters.N]
end


"""
    binary_tournament(population)
Apply binary tournament to obtain a solution from from population.
"""
function binary_tournament(population) 
    a, b = rand(1:length(population), 2)

    return is_better(population[a], population[b]) ? population[a] : population[b]
end


function binary_tournament(population, fitness::Vector{T})  where T <: Real
    a, b = rand(1:length(population), 2)

    return fitness[a] < fitness[b] ? population[a] : population[b]
end

mutable struct BiasedSelection
    num_elites::Int
    num_offsprings::Int
end


function selection(population, parameters::BiasedSelection)
    elites = 1:parameters.num_elites
    no_elites = parameters.num_elites+1:length(population)
    n = parameters.num_offsprings ÷ 2

    parent = ones(Int, 2n)
    parent[1:2:end] = rand(elites, n)
    parent[2:2:end] = rand(no_elites, n)
    parent
end
