##################################################################
# Source code taken from: https://github.com/jbytecode/MCCGA
# Credits to corresponding author.
##################################################################
include("utils.jl")

mutable struct MCCGA <: AbstractParameters
    N::Int # population size
    maxsamples::Int
    mutation::Float64
    probvector::Vector{Float64}
    use_local_search::Bool
end

"""
    MCCGA(;N, maxsamples)

### Parameters:

- `N` population size
- `maxsamples` maximum number of samples. 

"""
function MCCGA(;
        N = 100,
        maxsamples =10_000,
        mutation = 1 / N,
        use_local_search = true,
        information = Information(),
        options = Options()
    )

    parameters = MCCGA(N, maxsamples, mutation, [], use_local_search)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end

function initialize!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if options.iterations == 0
        options.iterations = 100_000_000
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end
    
    lower = problem.bounds[1,:]
    upper = problem.bounds[2,:]

    parameters.probvector = initialprobs(lower, upper, maxsamples = parameters.maxsamples)

    # sample a vector to create an initial State
    x = sample(parameters.probvector) |> floats
    initial_sol = create_solution(x, problem)
    return State(initial_sol, [initial_sol for i in 1:parameters.N])

end

function update_state!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    probvector = parameters.probvector
    chsize = length(probvector)

    # parents
    ch1 = sample(probvector)
    ch2 = sample(probvector)

    # evaluate cost function
    sol1 = create_solution(floats(ch1), problem)
    sol2 = create_solution(floats(ch2), problem)

    sol_winner = sol1
    winner = ch1
    loser  = ch2

    # check if ch2 (sol2) is better that ch1 (sol1)
    if is_better(sol2, sol1)
        sol_winner = sol2
        winner = ch2
        loser = ch1
    end

    # save in population the winner (informative only)
    # this is not used in the algorithm
    status.population[1 + status.iteration % parameters.N] = sol_winner

    # save best solution found so far
    is_better(sol_winner, status.best_sol) &&  (status.best_sol = sol_winner)

    for i = 1:chsize
        if winner[i] != loser[i]
            if winner[i] == 1
                probvector[i] += parameters.mutation
            else
                probvector[i] -= parameters.mutation
            end
            if probvector[i] > 1
                probvector[i] = 1
            elseif probvector[i] < 0
                probvector[i] = 0
            end
        end
    end
end


function final_stage!(
        status::State{xf_indiv}, # unconstrained case
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if !parameters.use_local_search
        status.final_time = time()
        return
    end

    costfunction(x) = evaluate(x, problem)

    probvector = parameters.probvector
    sampledvector = sample(probvector)

    # initial solution for the local search
    initial_solution = floats(sampledvector)

    options.debug && @info "Running NelderMead..."
    local_result = Optim.optimize(costfunction, initial_solution, Optim.NelderMead())
    # display Nelder-Mead result
    options.debug && display(local_result)
    options.debug && @info "NelderMead done!"

    # save best solution found so far!
    sol = create_child(local_result.minimizer, local_result.minimum)

    status.final_time = time()
    status.f_calls = problem.f_calls

    if is_better(sol, status.best_sol)
        status.best_sol = sol
    end

    return
end

function final_stage!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    status.final_time = time()
    # nothing to do for constrained or multi-objective
    @warn "MCCGA has been designed for unconstrained problems (local search not performed)."
end

function stop_criteria!(status, parameters::MCCGA, problem, information, options)
    # check budget limitation
    if status.stop
        return
    end

    mutation = parameters.mutation
    status.stop = all(x -> (x <= mutation) || (x >= 1.0 - mutation), parameters.probvector)

    status.stop && options.debug && @info "MCCGA stopping criteria is met."
end

