function final_stage!(
        status::State{xf_indiv}, # unconstrained case
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    status.final_time = time()
    !parameters.use_local_search &&  (return)
    if status.stop && status.termination_status_code ∈ [EVALUATIONS_LIMIT, TIME_LIMIT, ACCURACY_LIMIT]
        return
    end

    costfunction(x) = evaluate(x, problem)

    probvector = parameters.probvector
    sampledvector = sample(probvector)

    # initial solution for the local search
    initial_solution = floats(sampledvector)

    options.debug && @info "Running NelderMead..."
    local_result = Optim.optimize(costfunction, initial_solution, Optim.NelderMead())
    # display Nelder-Mead result
    options.debug && display(local_result)
    options.debug && @info "NelderMead done!"

    # save best solution found so far!
    sol = create_child(local_result.minimizer, local_result.minimum)

    status.final_time = time()
    status.f_calls = problem.f_calls

    if is_better(sol, status.best_sol)
        status.best_sol = sol
    end

    return
end
