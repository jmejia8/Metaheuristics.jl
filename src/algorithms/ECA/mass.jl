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

function getMass(U::Array{xf_indiv,1})
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness)
end

################################################
# Constrained
################################################
"""
    getMass(U, max_fitness, epsilon)
    get the mass in a vector `m[i] = f[i] + 2*max_fitness*sum_vio[i] `, tolerance
    in equality constraints are given by `epsilon`
"""
function getMass(U::Array{xfgh_indiv,1}, max_fitness, ε)
    n = length(U)
    fitness = zeros(Float64, n)

    for i = 1:n
        v = violationsSum(U[i].g, U[i].h, ε)
        if v > 0.0
            fitness[i] = U[i].f + 2*max_fitness*v
        else
            fitness[i] = U[i].f
        end
    end

    return fitnessToMass(fitness)

end


################################################
# Multi-objective
################################################
function getMass(U::Array{xFgh_indiv,1})
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        v = U[i].sum_violations
        if v > 0.0
            fitness[i] = v
        else
            fitness[i] = sum(U[i].f)
        end
    end

    return fitnessToMass(fitness)

end

