function default_stop_check(status, information, options)
    # budget termination
    call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)  ||
    time_stop_check(status, information, options) ||
    # accuracy
    accuracy_stop_check(status, information, options) ||
    # convergence
    stop_check(status, CheckConvergence(
                                        f_tol_abs = options.f_tol,
                                        f_tol_rel = options.f_tol_rel,
                                        x_tol =options.x_tol)) #=||
    stop_check(
               status,
               CheckConvergence(
                                criteria=[RobustConvergence(ftol=options.f_tol)]
                               )
              )
              =#
end
