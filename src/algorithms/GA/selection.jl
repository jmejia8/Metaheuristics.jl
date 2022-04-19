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

