function rank_solutions(solutions::Array{<:AbstractSolution}; is_better = is_better)
    mask = sortperm(solutions, lt=Metaheuristics.is_better, alg=Base.QuickSort)
    sortperm(mask, alg=Base.QuickSort)
end

function rank_solutions(solutions::Array{<: xFgh_solution}; is_better = is_better) 
    return non_dominated_sort(solutions)
end

get_fitness(population::Array{<:AbstractSolution}) = fvals(population)
get_fitness(population::Array{<: xFgh_solution}) = non_dominated_sort(solutions)

function _to_int_if_necessary(::Type{X}, v) where X<:Real
    return v
end

function _to_int_if_necessary(::Type{X}, v) where X<:Integer
    return round.(X, v)
end

