function _after_optimization!(status, parameters, problem, information, options, convergence)
    

    t = time()
    final_stage!(status, parameters, problem, information, options)
    status.overall_time += time() - t

    options.debug && @info "Termination reason: " * termination_status_message(status)
    # save convergence
    status.convergence = convergence;
end
