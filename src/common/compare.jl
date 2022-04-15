"""
    is_better(A, B)
    return true if A is better than B in a minimization problem.
    Feasibility rules and dominated criteria are used in comparison.
"""
function is_better(A::T, B::T) where T <: xf_solution
    return A.f < B.f
end

# feasibility rules
function is_better(A::T, B::T) where T <: xfgh_solution

    A_vio = A.sum_violations
    B_vio = B.sum_violations

    if A_vio < B_vio
        return true
    elseif B_vio < A_vio
        return false
    end


    return A.f < B.f

end

# is B dominated by A?
function is_better(A::T, B::T) where T <: xFgh_solution
    A_vio = A.sum_violations
    B_vio = B.sum_violations

    if A_vio < B_vio
        return true
    elseif B_vio < A_vio
        return false
    end

    
    compare(A.f, B.f) == 1

end


"""
    does A dominate B?
"""
dominates( A::T, B::T) where T <: xFgh_solution = is_better(A, B)


"""
    get_best(population)
    return best element in population according to the `is_better` function.
"""
function get_best(population::Array)
    if isempty(population)
        return nothing
    end
    
    best = population[1]
    for i in 2:length(population)
        if is_better(population[i], best)
            best = population[i]
        end
    end

    return best
end


"""
    compare(a, b)
    compares whether two vectors are dominated or not.
    Output:
    `1` if argument 1 (a) dominates argument 2 (b).
    `2` if argument 2 (b) dominates argument 1 (a).
    `3` if both arguments 1 (a) and 2 (b) are incomparable.
    `0` if both arguments 1 (a) and 2 (b) are equal.
"""
function compare(a::Vector, b::Vector)
    k = length(a)
    @assert k == length(b)

    i = 1
    while i <= k && a[i] == b[i]
        i += 1;
    end

    if i > k
        return 0 # equals
    end

    if a[i] < b[i]

        for j = i+1:k# (j = i+1; j <= k; ++j)
            if b[j] < a[j]
                return 3 #a and b are incomparable
            end
        end

        return 1; #  a dominates b
    end

    for j = i+1:k #(j = i+1; j < k; ++j) 
        if (a[j] < b[j])
            return 3 #; // a and b are incomparable
        end
    end

    return 2 # b dominates a
    
end


function compare(a::T, b::T) where T <: xFgh_solution

    A_vio = a.sum_violations
    B_vio = b.sum_violations

    if A_vio < B_vio
        return 1
    elseif B_vio < A_vio
        return 2
    end

    compare(a.f, b.f)
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

