abstract type AbstractVNS <: MH.AbstractParameters end


Base.@kwdef struct VND{I, N, L} <: AbstractVNS
    initial::I 
    neighborhood::N # neighborhood structures
    local_search::L  # local search strategy
end

function VND(;initial = nothing, neighborhood = nothing, local_search = nothing,
        options=MH.Options(), information=MH.Information())

    parameters = VND(initial, neighborhood, local_search)
    
    MH.Algorithm(parameters; options, information)
end


MH.iscompatible(::MH.BitArraySpace, ::AbstractVNS) = true
MH.iscompatible(::MH.PermutationSpace, ::AbstractVNS) = true

function MH.initialize!(status, parameters::AbstractVNS, problem, information, options, args...; kargs...)

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
        x = MH.minimizer(status)
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
        if MH.is_better(sol, status.best_sol)
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


function MH.final_stage!(status, parameters::AbstractVNS, problem, information, options, args...; kargs...)
end
