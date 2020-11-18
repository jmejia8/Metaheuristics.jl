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


function get_pareto_front(Population, is_better = is_better)
    ids = (Int[])
    n = length(Population)

    for i in 1:n
        s = true
        for j in 1:n
            i == j && continue

            if is_better(Population[j], Population[i])
                s = false
                break
            end
        end

        if s
            push!(ids, i)
        end
    end

    return Population[ids]

end

