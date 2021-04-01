################################################
# Unconstrained
################################################
function fitnessToMass(fitness::Vector{Float64})
    m = minimum(fitness)

    if m < 0
        fitness = 2 * abs(m) .+ fitness
    end
    
    fitness = 2 * maximum(fitness) .- fitness

    return fitness
end

"""
    getMass(U, max_fitness, epsilon)
    get the mass in a vector `m[i] = f[i] + 2*max_fitness*sum_vio[i] `, tolerance
    in equality constraints are given by `epsilon`
"""
function getMass(U::Array{xf_indiv,1})
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness)
end

