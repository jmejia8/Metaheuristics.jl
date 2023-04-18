"""
    gen_initial_state(problem,parameters,information,options)

Generate an initial state, i.e., compute uniformly distributed random vectors in bounds,
after that are evaluated in objective function. This method require that `parameters.N`
is valid attribute.
"""
gen_initial_state(problem,parameters,information,options,status::State{Any}) = gen_initial_state(problem,parameters,information,options)


function gen_initial_state(problem,parameters,information,options, status)
    if parameters.N != length(status.population)
        options.debug && @warn("Population size in provided State differs from that in parameters")
    end

    getdim(problem) != length(get_position(status.best_sol)) &&
        error("Invalid population (dimension does not match with bounds)")

    _complete_population!(status,problem,parameters,information,options)

    best_solution = get_best(status.population)
    # check if a better solution was found
    if is_better(status.best_sol, best_solution)
        best_solution = status.best_sol
    end

    return State(best_solution, status.population;f_calls = length(status.population))
end

function gen_initial_state(problem,parameters,information,options)
    # population array
    population = generate_population(parameters.N, problem, options.rng,ε=options.h_tol)

    # best solution
    best_solution = get_best(population)

    return State(best_solution, population; f_calls = length(population), iteration=1)
end


function _complete_population!(status,problem,parameters,information,options)
    if parameters.N < length(status.population)
        # increase population if necessary
        parameters.N = length(status.population)
        # TODO: use this options.debug == true to show the message?
        options.debug && @warn("Population size increased to $(parameters.N) due to initial solutions.")
        return
    end

    if parameters.N > length(status.population)
        # complete population with random in bounds
        n = parameters.N - length(status.population)
        missing_sols = generate_population(n, problem,ε=options.h_tol)
        # insert new solution into population
        append!(status.population, missing_sols)
    end
    
end

