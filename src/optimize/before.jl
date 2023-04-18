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
    status = method.status
    options.debug && @info("Initializing population...")
    start_time = time()
    # initialize population (... and optimizer)
    status = initialize!(status, parameters, problem, information, options)

    # save times
    status.start_time = start_time
    status.overall_time += time() - start_time
    method.status = status
    termination = method.termination

    # set termination based on convergence if necessary
    if isempty(termination.checkall) && isempty(termination.checkany)
        if fval(first(status.population)) isa Array
            # multi-objective termination criterion
            c = RobustConvergence(ftol=options.f_tol)
            options.debug && @info("Set termination criteria: RobustConvergence")
        else
            # single-objective criteria (convergence based)
            c = CheckConvergence(f_tol_abs = options.f_tol,
                                 f_tol_rel = options.f_tol_rel,
                                 x_tol = options.x_tol)
            options.debug && @info("Set termination criteria: convergence indicators")
        end
        
        push!(termination.checkany, c)
    end
    
    
    # show status if necessary
    show_status(status, parameters, options)

    status.iteration = 1
    convergence = State{typeof(status.best_sol)}[]
    options.store_convergence && update_convergence!(convergence, status)

    options.debug && @info("Starting main loop...")
    logger(status)

    status.stop = status.stop ||
                    # common stop criteria (budget based)
                    default_stop_check(status, information, options) ||
                    # stop criteria given by user
                    stop_check(status, termination)

    status, parameters, problem, information, options, convergence
end

