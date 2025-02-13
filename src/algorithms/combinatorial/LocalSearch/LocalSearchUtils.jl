abstract type AbstractLocalSearch end
struct BestImproveSearch  <: AbstractLocalSearch end
struct FirstImproveSearch <: AbstractLocalSearch end

include("neighborhood.jl")

function local_search(x, neighbourhood::Neighborhood, ls::AbstractLocalSearch, problem)
    # creates an iterator and perform local search over the iterator
    iter = NeighborhoodIterator(x, neighbourhood)
    local_search(x, iter, ls, problem)
end

function local_search(x, neighbourhood::InternalNeighborhood, ::BestImproveSearch, problem)
    best = create_solution(copy(x), problem)
    for xnew in neighbourhood
        sol = create_solution(xnew, problem)
        if is_better(sol, best)
            best = deepcopy(sol)
        end
    end
    best
end

function local_search(x, neighbourhood::InternalNeighborhood, ::FirstImproveSearch, problem)
    initial = create_solution(copy(x), problem) 
    for xnew in neighbourhood
        sol = create_solution(xnew, problem)
        if is_better(sol, initial)
            return deepcopy(sol)
        end
    end
    initial
end

function local_search(x, ls::AbstractLocalSearch, problem)
    if problem.search_space isa PermutationSpace
        neighbourhood = TwoOptNeighborhood()
    elseif problem.search_space isa BoxConstrained
        neighbourhood = GridSampler(problem.search_space)
    elseif problem.search_space isa BitArraySpace
        # TODO
        neighbourhood = nothing
    else
        @error "Definition of local_search for $(problem.search_space) is required."
        return
    end
    local_search(x, neighbourhood, ls, problem)
end


