abstract type AbstractGRASP <: MH.AbstractParameters end

struct GRASP{I, T, L} <: AbstractGRASP
	initial::I
	constructor::T
	local_search::L
end

Base.@kwdef struct GreedyRandomizedContructor
    candidates
    instance = nothing
    α::Float64 = 0.4
    rng = MH.default_rng_mh()
end


function construct(constructor)
    @error "Define your own randomized greedy constructor"
    status.stop = true
    nothing
end

function local_search(x, localsearch, problem)
    @error "Define your own local search"
    nothing
end


function compute_cost(candidates, constructor, instance)
    @warn "Define compute_cost for\nconstructor=$constructor\ninstance=$instance"
    zeros(length(candidates))
end


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
        # select candidate at random and insert into solution
        s = rand(constructor.rng, RCL)
        push!(S, candidates[s])
        # update list of candidates
        deleteat!(candidates, s)
    end
    S
end

function GRASP(;initial=nothing, constructor=nothing, local_search=nothing,
                options = MH.Options(), information=MH.Information())
	# TODO
	if isnothing(constructor) || isnothing(local_search)
        error("Provide a constructor and a local search")
	end
    grasp = GRASP(initial, constructor, local_search)
    MH.Algorithm(grasp; options, information)
end

MH.iscompatible(::MH.BitArraySpace, ::AbstractGRASP) = true
MH.iscompatible(::MH.PermutationSpace, ::AbstractGRASP) = true

function MH.initialize!(status, parameters::AbstractGRASP, problem, information, options, args...; kargs...)

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
    
    
    sol = MH.create_solution(x0, problem)
	# TODO
	MH.State(sol, [sol])
end

function MH.update_state!(
        status,
        parameters::AbstractGRASP,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    # heuristic construction
    x = construct(parameters.constructor)
    # perform local search and evaluate solutions
    x_improved = local_search(x, parameters.local_search, problem)

    # since local search can return a solution without evaluation
    # it is necessary to evaluate objective function
    if x_improved isa AbstractVector
        # evaluate solution
        sol = MH.create_solution(x_improved, problem)
    elseif x_improved isa MH.AbstractSolution
        sol = x_improved
    else
        # seems that local search returned something different to a vector
        return
    end
    
    # save best solutions
    if MH.is_better(sol, status.best_sol)
        status.best_sol = sol
    end

    # update history (for convergence checking)
    #= TODO
    push!(status.population, sol)
    max_history = 10
    if length(status.population) > 2max_history
        offprings = status.population[max_history:end]
        deleteat!(status.population, max_history:length(status.population))
        MH.environmental_selection!(status.population, offprings, MH.ElitistReplacement())
    end
    =#

end

function MH.final_stage!(
        status,
        parameters::AbstractGRASP,
        problem,
        information,
        options,
        args...;
        kargs...
   )
   # TODO
   nothing
end

