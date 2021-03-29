function initialize!(problem,engine,parameters,status,information,options)
    a, b = problem.bounds[1,:], problem.bounds[2,:]
    D = length(a)

    # population array
    Population = generate_population(problem.f, parameters.N, problem.bounds,ε=options.h_tol)
    status.population = Population
    # current evaluations
    status.f_calls = parameters.N


    # current generation
    status.iteration = 0

    # best solution
    status.best_sol = get_best(Population)

    # status.stop = engine.stop_criteria(status, information, options)

end

function gen_initial_state(problem,parameters,information,options)
    # population array
    population = generate_population(problem.f, parameters.N, problem.bounds,ε=options.h_tol)

    # best solution
    best_solution = get_best(population)

    return State(best_solution, population; f_calls = length(population), iteration=1)
end

