function fast_non_dominated_sort!(pop, is_better)
    n = length(pop)

    dom_list = [ Int[] for i in 1:n ]
    rank = zeros(Int, n)
    dom_count = zeros(Int, n)

    for i in 1:n
        for j in i+1:n
            comparison = compare(pop[i].f, pop[j].f)
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

    for i in 1:length(pop)
        pop[i].rank = rank[i]
    end

    sort!(pop, by = s -> s.rank, alg = Base.Sort.QuickSort)
    pop
end



function update_crowding_distance!(pop)

    for sol in pop
        sol.crowding = 0.0
    end


    if length(pop[1].f) == 2
        sort!(pop, by = x -> x.f[1])
        pop[1].f[1] == pop[end].f[1] && return #Don't waste time if all indivs are the same
        pop[1].crowding = pop[end].crowding = Inf

        width_y1 = (pop[end].f[1] - pop[1].f[1])
        width_y2 = (pop[1].f[2] - pop[end].f[2])
        for i = 2:length(pop)-1
            pop[i].crowding = (pop[i+1].f[1] - pop[i-1].f[1]) / width_y1 + (pop[i-1].f[2] - pop[i+1].f[2]) / width_y2
        end

        return pop
    end

    for j = 1:length(first(pop).f) # Foreach objective
        let j = j
            sort!(pop, by = x -> x.f[j]) #sort by the objective value
        end
        pop[1].crowding = pop[end].crowding = Inf #Assign infinite value to extremas
        if pop[1].f[j] != pop[end].f[j]
            for i = 2:length(pop)-1
                pop[i].crowding += (pop[i+1].f[j] - pop[i-1].f[j]) / (pop[end].f[j] - pop[1].f[j])
            end
        end
    end
    pop
end



function truncate_population!(population, N, is_better)
    fast_non_dominated_sort!(population, is_better)

    #update_crowding_distance!(population)

    let f::Int = 1
        ind = 0
        indnext = findlast(x -> x.rank == f, population)
        while 0 < indnext <= N
            ind = indnext
            f += 1
            indnext = findlast(x -> x.rank == f, population)
        end
        indnext == 0 && (indnext = length(population))
        update_crowding_distance!(view(population, ind+1:indnext))
        sort!(view(population, ind+1:indnext), by = x -> x.crowding, rev = true, alg = PartialQuickSort(N-ind))
    end
    deleteat!(population, N+1:length(population))
end
