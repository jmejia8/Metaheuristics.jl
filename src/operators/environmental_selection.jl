"""
    GenerationalReplacement()

Generational replacement.
"""
struct GenerationalReplacement end

function environmental_selection!(population, offsprings, parameters::GenerationalReplacement)
    @assert length(population) == length(offsprings)
    population[:] = offsprings
end

"""
    ElitistReplacement()

Offspring is inserted in population to keep the best individuals (keep population size).
"""
struct ElitistReplacement end

function environmental_selection!(
        population,
        offsprings,
        parameters::ElitistReplacement;
        is_better = is_better
    )

    N = length(population)
    append!(population, offsprings)
    sort!(population, lt=is_better)
    deleteat!(population, N+1:length(population))
end


"""
    RankAndCrowding()

Perform `environmental_selection!` based non-dominated ranking and crowding
distance.
"""
struct RankAndCrowding end

function environmental_selection!(
        population,
        offsprings,
        parameters::RankAndCrowding;
        is_better = is_better
    )

    N = length(population)
    append!(population, offsprings)
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
