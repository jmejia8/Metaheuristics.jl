function compute_crowding_distance(pop)
    crowding = zeros(length(pop))


    if length(pop[1].f) == 2
        sort!(pop, by = x -> x.f[1])
        pop[1].f[1] == pop[end].f[1] && (return crowding) #Don't waste time if all indivs are the same
        crowding[1] = crowding[end] = Inf

        width_y1 = (pop[end].f[1] - pop[1].f[1])
        width_y2 = (pop[1].f[2] - pop[end].f[2])
        for i = 2:length(pop)-1
            crowding[i] = (pop[i+1].f[1] - pop[i-1].f[1]) / width_y1 + (pop[i-1].f[2] - pop[i+1].f[2]) / width_y2
        end

        return crowding
    end

    for j = 1:length(first(pop).f) # Foreach objective
        let j = j
            sort!(pop, by = x -> x.f[j]) #sort by the objective value
        end
        crowding[1] = crowding[end] = Inf #Assign infinite value to extremas
        if pop[1].f[j] != pop[end].f[j]
            for i = 2:length(pop)-1
                crowding[i] += (pop[i+1].f[j] - pop[i-1].f[j]) / (pop[end].f[j] - pop[1].f[j])
            end
        end
    end

    crowding
end

function update_crowding_distance!(pop)
    # at leat 3 items to compute crowding distance
    length(pop) < 3 && (return pop) 
    
    crowding = compute_crowding_distance(pop)
    for (i, sol) in enumerate(pop)
        sol.crowding = crowding[i]
    end

    pop
end



function truncate_population!(population, N)
    fast_non_dominated_sort!(population)

    let f::Int = 1
        ind = 0
        indnext = findlast(x -> x.rank == f, population)
        while !isnothing(indnext) && 0 < indnext <= N
            ind = indnext
            f += 1
            indnext = findlast(x -> x.rank == f, population)
        end
        isnothing(indnext) && (indnext = length(population)) 
        update_crowding_distance!(view(population, ind+1:indnext))
        sort!(view(population, ind+1:indnext), by = x -> x.crowding, rev = true, alg = PartialQuickSort(N-ind))
    end
    deleteat!(population, N+1:length(population))
end
