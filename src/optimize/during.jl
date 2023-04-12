function _during_optimization!(status, parameters, problem, information, options, convergence, logger, termination)
    status.iteration += 1
    t = time()
    # perform one optimization step
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
    status.stop = status.stop ||
                    # common stop criteria (budget based)
                    default_stop_check(status, information, options) ||
                    # stop criteria given by user
                    stop_check(status, termination)

    # TODO: check if following method should be deprecated in v4
    # stop criteria defined by optimizer
    !status.stop && stop_criteria!(status, parameters, problem, information, options)
end
