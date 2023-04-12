function _during_optimization!(status, parameters, problem, information, options, convergence, logger)
    status.iteration += 1
    # perform one optimization step
    t = time()
    update_state!(status, parameters, problem, information, options)
    status.overall_time += time() - t
    # store the number of function evaluations
    status.f_calls = problem.f_calls
    # show status if debug = true
    show_status(status, parameters, options)
    # store convergence if necessary
    options.store_convergence && update_convergence!(convergence, status)
    # user defined logger
    logger(status)
    # common stop criteria
    status.stop = status.stop || default_stop_check(status, information, options)

    # user defined stop criteria
    stop_criteria!(status, parameters, problem, information, options)
end
