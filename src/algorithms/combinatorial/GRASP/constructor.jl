"""
    GreedyRandomizedContructor(;candidates, instance, α, rng)

This structure can be used to use the Greedy Randomized Contructor.

- `candidates` vector of candidates (item to choose to form a solution).
- `instance` a user define structure for saving data on the problem instance.
- `α`  controls the randomness (0 deterministic greedy heuristic, 1 for pure random).
- `rng` storages the random number generator.

This constructor assumes overwriting the `compute_cost` function:

```julia
import Metaheuristics as MH
struct MyInstance end
function MH.compute_cost(candidates, constructor, instance::MyInstance)
    # ...
end
```

See also [`compute_cost`](@ref) and [`GRASP`](@ref)
"""
Base.@kwdef struct GreedyRandomizedContructor
    candidates
    instance = nothing
    α::Float64 = 0.6
    rng = default_rng_mh()
end

"""
    compute_cost(candidates, constructor, instance)

Compute the cost for each candidate in `candidates`, for given `constructor` and
provided `instance`.

See also [`GreedyRandomizedContructor`](@ref) and [`GRASP`](@ref)
"""
function compute_cost(candidates, constructor, instance)
    @warn "Define compute_cost for\nconstructor=$constructor\ninstance=$instance"
    zeros(length(candidates))
end

"""
    construct(constructor)

Constructor procedure for GRASP.

See also [`GreedyRandomizedContructor`](@ref), [`compute_cost`](@ref) and [`GRASP`](@ref)
"""
function construct(constructor::GreedyRandomizedContructor)
    candidates = constructor.candidates |> copy
    α = constructor.α
    # create empty solution S
    S = empty(candidates)
    # construct solution
    while !isempty(candidates)
        cost = compute_cost(candidates, constructor, constructor.instance)
        cmin = minimum(cost)
        cmax = maximum(cost)
        # compute restricted candidate list
        RCL = [i for i in eachindex(candidates) if cost[i] <= cmin + α*(cmax - cmin) ]
        if isempty(RCL)
            @error "RCL is empty. Try increasing α or check your `compute_cost` method."
            return
        end
        
        # select candidate at random and insert into solution
        s = rand(constructor.rng, RCL)
        push!(S, candidates[s])
        # update list of candidates
        deleteat!(candidates, s)
    end
    S
end
