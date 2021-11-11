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
