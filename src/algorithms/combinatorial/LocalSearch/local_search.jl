abstract type AbstractLocalSearch end
struct BestImprovingSearch  <: AbstractLocalSearch end
struct FirstImprovingSearch <: AbstractLocalSearch end

include("neighbourhood.jl")

function local_search(x, neighbourhood::Neighbourhood, ls::AbstractLocalSearch, problem)
    # creates an iterator and perform local search over the iterator
    iter = NeighborhoodIterator(x, neighbourhood)
    local_search(x, iter, ls, problem)
end

function local_search(x, neighbourhood::InternalNeighbourhood, ::BestImprovingSearch, problem)
    best = MH.create_solution(copy(x), problem)
    for xnew in neighbourhood
        sol = MH.create_solution(xnew, problem)
        if MH.is_better(sol, best)
            best = deepcopy(sol)
        end
    end
    best
end

function local_search(x, neighbourhood::InternalNeighbourhood, ::FirstImprovingSearch, problem)
    initial = MH.create_solution(copy(x), problem) 
    for xnew in neighbourhood
        sol = MH.create_solution(xnew, problem)
        if MH.is_better(sol, initial)
            return deepcopy(sol)
        end
    end
    initial
end

function local_search(x, ls::AbstractLocalSearch, problem)
    if problem.search_space isa MH.PermutationSpace
        neighbourhood = TwoOpt(;x)
    elseif problem.search_space isa MH.BoxConstrained
        neighbourhood = MH.GridSampler(problem.search_space)
    elseif problem.search_space isa MH.BitArraySpace
        # TODO
        neighbourhood = nothing
    else
        @error "Definition of local_search for $(problem.search_space) is required."
        return
    end
    local_search(x, neighbourhood, ls, problem)
end


