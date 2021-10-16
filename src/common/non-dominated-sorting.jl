"""
    non_dominated_sort(population)

Return the non-dominated rank from population.
"""
function non_dominated_sort(pop)
    n = length(pop)

    dom_list = [ Int[] for i in 1:n ]
    rank = zeros(Int, n)
    dom_count = zeros(Int, n)

    for i in 1:n
        for j in i+1:n
            comparison = compare(pop[i], pop[j])
            if comparison == 1 #is_better(pop[i], pop[j])
                push!(dom_list[i], j)
                dom_count[j] += 1
            elseif comparison == 2 #is_better(pop[j], pop[i])
                push!(dom_list[j], i)
                dom_count[i] += 1
            end
        end
        if dom_count[i] == 0
            rank[i] = 1
        end
    end

    k = UInt16(2)
    while any(==(k-one(UInt16)), (rank[p] for p in 1:n)) #ugly workaround for #15276
        for p in 1:n
            if rank[p] == k-one(UInt16)
                for q in dom_list[p]
                    dom_count[q] -= one(UInt16)
                    if dom_count[q] == zero(UInt16)
                        rank[q] = k
                    end
                end
            end
        end
        k += one(UInt16)
    end

    return rank
end

"""
    get_fronts(population, computed_ranks = true)

Return each sub-front in an array. If `computed_ranks == true`, this method assumes that
`fast_non_dominated_sort!(population)` has been called before.
"""
function get_fronts(population, computed_ranks = true)
    if computed_ranks
        ranks = [s.rank for s in population]
    else
        ranks = non_dominated_sort(population)
    end

    fronts = [Int[] for i in eachindex(unique(ranks))]

    for (i, r) in enumerate(ranks)
        push!(fronts[r], i)
    end

    return fronts
end


function fast_non_dominated_sort!(pop)
    rank = non_dominated_sort(pop)

    for i in 1:length(pop)
        pop[i].rank = rank[i]
    end

    sort!(pop, by = s -> s.rank, alg = Base.Sort.QuickSort)
    pop
end


function get_non_dominated_solutions_perm(population)
    ids = Int[1]
    n = length(population)

    for i in 2:n
        j = 1
        while j <= length(ids)
            jj = ids[j]
            relation = compare(population[i], population[jj])
 
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

    return ids
end


function get_non_dominated_solutions(population)
    return population[get_non_dominated_solutions_perm(population)]
end
