struct UniformCrossover
    p::Float64
    UniformCrossover(;p=0.5) = new(p)
end

function crossover(population, parameters::UniformCrossover)
    n = length(population) ÷ 2
    offspring_A = positions(population[1:n])
    offspring_B = positions(population[n+1:2n])
    mask = rand(size(offspring_A)...) .<= parameters.p
    tmp = copy(offspring_A[mask])
    offspring_A[mask] = offspring_B[mask]
    offspring_B[mask] = tmp
    [offspring_A; offspring_B]
end


struct OrderCrossover end

function crossover(population, parameters::OrderCrossover)
    O = positions(population)
    N, D = size(O)
    s = rand(1:D, N) # slash points
    for i = 1:2:N
        PA = O[i, :]   # parent A
        PB = O[i+1, :] # parent B
        O[i,  s[i]+1:D] = setdiff(PB, PA[1:s[i]]);
        O[i+1,s[i]+1:D] = setdiff(PA, PB[1:s[i]]);
    end
    O
end
