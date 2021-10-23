"""
    binary_tournament(population)
Apply binary tournament to obtain a solution from from population.
"""
function binary_tournament(population) 
    a, b = rand(1:length(P), 2)

    return is_better(population[a], population[b]) ? population[a] : population[b]
end
