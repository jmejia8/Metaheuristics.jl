function initialize!(problem,engine,parameters,status,information,options)
    a, b = problem.bounds[1,:], problem.bounds[2,:]
    D = length(a)

    # population array
    Population = initializePop(problem.f, parameters.N, D, a, b)
    status.population = Population
    # current evaluations
    status.f_calls = parameters.N


    # current generation
    status.iteration = 0

    # best solution
    status.best_sol = getBest(Population, :minimize)

    status.stop = engine.stop_criteria(status, information, options)
    
end