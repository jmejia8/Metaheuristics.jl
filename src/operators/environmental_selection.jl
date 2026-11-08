"""
    GenerationalReplacement()

Generational replacement.
"""
struct GenerationalReplacement end

function environmental_selection!(
    population,
    offsprings,
    parameters::GenerationalReplacement,
)
    if length(population) != length(offsprings)
        throw(DimensionMismatch("population and offsprings must have the same length."))
    end
    copyto!(population, offsprings)
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
    is_better = is_better,
)

    N = length(population)
    append!(population, offsprings)
    sort!(population, lt = is_better)
    deleteat!(population, N + 1:length(population))
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
    is_better = is_better,
)

    N = length(population)
    append!(population, offsprings)
    truncate_population!(population, N)
end

"""
    truncate_population!(population, N)

Keep the best `N` solutions in `population` according to non-dominated ranking and
crowding distance (as in NSGA-II).
"""
function truncate_population!(population, N)
    fast_non_dominated_sort!(population)
    L = length(population)
    L == 0 && return population

    # population is sorted by rank: find the last front that fits in N (ends at ind)
    # and the end of the next front (indnext)
    ind = 0
    indnext = _front_end(population, 1)
    while indnext <= N && indnext < L
        ind = indnext
        indnext = _front_end(population, ind + 1)
    end
    indnext <= N && (ind = indnext) # every solution fits

    update_crowding_distance!(view(population, ind + 1:indnext))
    sort!(
        view(population, ind + 1:indnext),
        by = x -> x.crowding,
        rev = true,
        alg = PartialQuickSort(N - ind),
    )
    deleteat!(population, N + 1:length(population))
end

# last index of the front that starts at `start` (population sorted by rank)
function _front_end(population, start)
    r = population[start].rank
    e = start
    @inbounds while e < length(population) && population[e + 1].rank == r
        e += 1
    end
    e
end
