abstract type AbstractVNS <: AbstractParameters end


struct VND{I, N, L} <: AbstractVNS
    initial::I 
    neighborhood::N # neighborhood structures
    local_search::L  # local search strategy
end

"""
    VND(;initial, neighborhood, local_search,...)

Variable Neighborhood Descent.

# Allowed parameters

- `initial`: Use this parameter to provide an initial solution (optional).
- `neighborhood`: Neighborhood structure.
- `local_search` the local search strategy `BestImproveSearch()` (default) and `FirstImproveSearch()`.


# Example: Knapsack Problem

```julia
import Metaheuristics as MH

struct MyKPNeighborhood <: MH.Neighborhood
    k::Int
end

function MH.neighborhood_structure(x, s::MyKPNeighborhood, i::Integer)
    # return the i-th neighbor around x, regarding s.k structure
    i > length(x) && return nothing
    reverse!(view(x, i:min(length(x), i+s.k)))
    x
end


function main()
    profit = [55, 10,47, 5, 4, 50, 8, 61, 85, 87]
    weight = [95, 4, 60, 32, 23, 72, 80, 62, 65, 46]
    capacity = 269

    # objective function and search space
    f, search_space, _ = MH.TestProblems.knapsack(profit, weight, capacity)

    # list the neighborhood structures
    neighborhood = [MyKPNeighborhood(1), MyKPNeighborhood(2), MyKPNeighborhood(3)]
    local_search = MH.BestImproveSearch()
    # instantiate VNS
    vnd = MH.VND(;neighborhood, local_search)

    res = MH.optimize(f, search_space, vnd)
    display(res)
end

main()
```
"""
function VND(;initial = nothing, neighborhood = nothing, local_search = BestImproveSearch(),
        options=Options(), information=Information())

    parameters = VND(initial, neighborhood, local_search)
    
    Algorithm(parameters; options, information)
end


iscompatible(::BitArraySpace, ::AbstractVNS) = true
iscompatible(::PermutationSpace, ::AbstractVNS) = true

function initialize!(status, parameters::AbstractVNS, problem, information, options, args...; kargs...)

    if isnothing(parameters.initial)
        x0 = rand(options.rng, problem.search_space)
    else
        x0 = parameters.initial
    end
    # set default budget
    options.f_calls_limit = Inf
    if options.iterations <= 0
        options.iterations = 500
    end
    
    sol = create_solution(x0, problem)
	# TODO
	State(sol, [sol])
end

function update_state!(
        status,
        parameters::VND,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    improvement = false 
    l = 1

    #  check if movement is required
    while l <= length(parameters.neighborhood)
        # current solution
        x = minimizer(status)
        # exploration of the neighborhood
        neighborhood = parameters.neighborhood[l]
        # local search around x
        sol = local_search(x, neighborhood, parameters.local_search, problem)

        # check for empty neighborhood
        if isnothing(sol)
            l += 1
            continue
        end

        # move or not
        if is_better(sol, status.best_sol)
            status.best_sol = sol
            l = 1
            improvement = true
        else
            l += 1
        end
    end

    # stop if no improvement is obtained
    status.stop = !improvement
end


function final_stage!(status, parameters::AbstractVNS, problem, information, options, args...; kargs...)
end
