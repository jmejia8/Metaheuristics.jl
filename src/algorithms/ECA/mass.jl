function fitnessToMass(fitness::Vector{Float64}, searchType::Symbol)
    m = minimum(fitness)

    if m < 0
        fitness = 2 * abs(m) .+ fitness
    end

    if searchType == :minimize
        fitness = 2 * maximum(fitness) .- fitness
    end

    return fitness
end

function getMass(U::Array{xf_indiv,1}, searchType::Symbol; kargs...)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness, searchType)
end

function getMass(U::Array{xfg_indiv,1}, searchType; kargs...)
    return getMass(U, searchType)

end

function getMass(U::Array{xfgh_indiv,1}, searchType; ε = 0.0)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        v = U[i].sum_violations
        if v > 0.0
            fitness[i] = v
        else
            fitness[i] = U[i].f
        end
    end

    return fitnessToMass(fitness, searchType)

end

function getMass(U::Array{xFgh_indiv,1}, searchType; ε = 0.0)
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

    return fitnessToMass(fitness, searchType)

end

