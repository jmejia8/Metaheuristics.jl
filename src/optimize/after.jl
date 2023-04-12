function _after_optimization!(status, parameters, problem, information, options, convergence)
    
    status.overall_time = time() - status.start_time

    final_stage!(status, parameters, problem, information, options)

    options.debug && @info "Termination reason: " * termination_status_message(status)
    # save convergence
    status.convergence = convergence;
end
