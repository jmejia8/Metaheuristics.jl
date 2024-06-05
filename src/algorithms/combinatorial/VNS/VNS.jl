include("VND.jl")

struct VNS{I, S, L, C, N} <: AbstractVNS
    initial::I
    neighborhood_shaking::S
    neighborhood_local::L
    local_search::C
    neighborhood_change::N
end

# change neighborhood
struct SequentialChange end
struct CyclicChange end

"""
    VNS(;initial, neighborhood_shaking, neighborhood_local, local_search, neighborhood_change)

General Variational Neighborhood Search.

# Allowed parameters

- `initial`: Use this parameter to provide an initial solution (optional).
- `neighborhood_shaking`: Neighborhood structure for the shaking step.
- `neighborhood_local`: Neighborhood structure for the local search.
- `local_search`: the local search strategy `BestImproveSearch()` (default) and `FirstImproveSearch()`.
- `neighborhood_change`: The procedure for changing among neighborhood structures  (default `SequentialChange()`).


# Example: Knapsack Problem

```julia
import Metaheuristics as MH

struct MyKPNeighborhood <: MH.Neighborhood
    k::Int
end

function MH.neighborhood_structure(x, s::MyKPNeighborhood, rng)
    # this is defined due to shaking procedure requires a random one
    # not the i-th neighbor.
    i = rand(rng, 1:length(x))
    reverse!(view(x, i:min(length(x), i+s.k)))
    x
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
    nitems = length(weight)

    # objective function and search space
    f, search_space, _ = MH.TestProblems.knapsack(profit, weight, capacity)

    # list the neighborhood structures
    neighborhood_shaking = [MyKPNeighborhood(6), MyKPNeighborhood(5), MyKPNeighborhood(4)]
    neighborhood_local = [MyKPNeighborhood(3), MyKPNeighborhood(2), MyKPNeighborhood(1)]
    local_search = MH.BestImproveSearch()
    # instantiate VNS
    vnd = MH.VNS(;neighborhood_shaking, neighborhood_local, local_search, options=MH.Options(verbose=true))

    res = MH.optimize(f, search_space, vnd)
    display(res)
end

main()
```
"""
function VNS(;initial=nothing,neighborhood_shaking=nothing, neighborhood_local=nothing,
        local_search=FirstImproveSearch(), neighborhood_change=SequentialChange(),
        options=Options(), information=Information())

    # TODO
    if isnothing(neighborhood_shaking) && isnothing(neighborhood_local) 
        error("Provide neighborhood_shaking and neighborhood_local.")
    end

    parameters = VNS(initial, neighborhood_shaking, neighborhood_local,
                     local_search, neighborhood_change)

    Algorithm(parameters; options, information)
end

function shake(x, neighborhood, rng)
    # select at random
    neighborhood_structure(x, neighborhood, rng)
end

function neighborhood_change(old, new, k, ::SequentialChange)
    if is_better(new, old)
        return 1, new
    end
    k + 1, old
end

function neighborhood_change(old, new, k, ::CyclicChange)
    k += 1
    if is_better(new, old)
        return k, new
    end
    k, old
end

function update_state!(status, parameters::VNS, problem, information, options, args...; kargs...)
    # current solution
    sol = first(status.population)
    k = 1
    while k <= length(parameters.neighborhood_shaking)
        x = get_position(sol)
        neighborhood = parameters.neighborhood_shaking[k]
        # select x at random from kth neighborhood
        xp = shake(x, neighborhood, options.rng)

        # perform local search around xp using VND
        # TODO: update this for considering other VNS variants (for the local search)
        vnd = VND(;initial=xp, neighborhood=parameters.neighborhood_local,
                  local_search = parameters.local_search)
        _res_local = optimize(problem.f, problem.search_space, vnd)
        sol_new = _res_local.best_sol

        # neighborhood change or not?
        k, sol = neighborhood_change(sol, sol_new, k, parameters.neighborhood_change)

        # save best result so far (internal use only, VNS doesn't use it)
        if is_better(sol_new, status.best_sol)
            status.best_sol = sol_new
        end
        problem.f_calls += _res_local.f_calls

    end
    # save sol (x) for the next iteration of VNS
    status.population = [sol]
end

