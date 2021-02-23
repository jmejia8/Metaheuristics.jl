"""
    is_better(A, B)
    return true if A is better than B in a minimization problem.
    Feasibility rules and dominated criteria area used in comparison.
"""
function is_better(A::xf_indiv, B::xf_indiv)
    return A.f < B.f
end

# feasibility rules
function is_better(A::xfgh_indiv, B::xfgh_indiv)

    A_vio = A.sum_violations
    B_vio = B.sum_violations

    if A_vio < B_vio
        return true
    elseif B_vio > A_vio
        return false
    end


    return A.f < B.f

end

# is B dominated by A?
function is_better(
    A::xFgh_indiv,
    B::xFgh_indiv
)

    A_vio = A.sum_violations
    B_vio = B.sum_violations

    if A_vio < B_vio
        return true
    elseif B_vio > A_vio
        return false
    end

    
    # A dominates an objective in B
    A_dominates_B = false
    for i in 1:length(A.f)
        if B.f[i] < A.f[i]
            return false
        elseif A.f[i] < B.f[i]
            A_dominates_B = true
        end
    end

    return A_dominates_B

end


"""
    does A dominate B?
"""
dominates( A::xFgh_indiv, B::xFgh_indiv) = is_better(A, B)


"""
    get_best(population)
    return best element in population according to the `is_better` function.
"""
function get_best(population::Array)
    if isempty(population)
        return nothing
    end
    
    best = population[1]
    for sol in population
        if is_better(sol, best)
            best = sol
        end
    end

    return best
end

"""
    argworst(population)
    return the index of the worst element in population
"""
function argworst(population::Array)
    if isempty(population)
        return -1
    end
    
    i_worst = 1
    for i in 2:length(population)
        if is_better(population[i_worst], population[i])
            i_worst = i
        end
    end

    return i_worst
end


"""
    argworst(population)
    return the index of the worst element in population
"""
function argbest(population::Array)
    if isempty(population)
        return -1
    end
    
    i_best = 1
    for i in 2:length(population)
        if is_better(population[i], population[i_best])
            i_best = i
        end
    end

    return i_best
end

