function get_non_dominated_solutions(Population, is_better = is_better)
    ids = Int[1]
    n = length(Population)

    for i in 2:n
        j = 1
        while j <= length(ids)
            jj = ids[j]
            relation = compare(Population[i], Population[jj])
 
            if relation == 2 # j dominates i
                break
            elseif relation == 1 # i dominates j
                deleteat!(ids, j)
                continue
            end

            j += 1
        end

        if j > length(ids)
            push!(ids, i)
        end
        
    end

    return Population[ids]

end

