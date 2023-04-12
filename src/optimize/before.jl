function _before_optimization!(f, search_space, method, logger)
    # check whether optimizer is compatible with search_space
    check_compat(search_space, method.parameters)

    # TODO: method = deepcopy(method)
    information = method.information
    options = method.options
    parameters = method.parameters
    problem = Problem(f, search_space; parallel_evaluation=options.parallel_evaluation)
    seed!(options.seed)

    # initialization steps
    start_time = time()
    status = method.status
    options.debug && @info("Initializing population...")
    status = initialize!(status,parameters, problem, information, options)
    status.start_time = start_time
    method.status = status
    
    # show status if necessary
    show_status(status, parameters, options)

    status.iteration = 1
    convergence = State{typeof(status.best_sol)}[]
    options.store_convergence && update_convergence!(convergence, status)

    options.debug && @info("Starting main loop...")
    logger(status)
    status.stop = status.stop || default_stop_check(status, information, options)

    status, parameters, problem, information, options, convergence
end

