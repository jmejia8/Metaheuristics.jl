"""
    call_limit_stop_check(status, information, options)

Limit the number of function evaluations, i.e., 
return `status.f_calls >= options.f_calls_limit`.
"""
function call_limit_stop_check(status, information, options)
    cond =  status.f_calls >= options.f_calls_limit
    cond && (status.termination_status_code = EVALUATIONS_LIMIT)
    cond
end

"""
    iteration_stop_check(status, information, options)
Used to limit the number of iterations.
"""
function iteration_stop_check(status, information, options)
    cond =  status.iteration >= options.iterations
    cond && (status.termination_status_code = ITERATION_LIMIT)
    cond
end

"""
    time_stop_check(status, information, options)
Used to limit the time (in seconds), i.e., `status.overall_time >= options.time_limit`.
"""
function time_stop_check(status, information, options)
    cond =  status.overall_time >= options.time_limit
    cond && (status.termination_status_code = TIME_LIMIT)
    cond
end

"""
    accuracy_stop_check(status, information, options)

If the optimum is provided, then check if the accuracy is met via
`abs(status.best_sol.f - information.f_optimum) < options.f_tol`.
"""
function accuracy_stop_check(status::State{<:xf_solution}, information, options)
    cond =  !isnan(information.f_optimum) && abs(fval(status.best_sol) - information.f_optimum) < options.f_tol
    cond && (status.termination_status_code = ACCURACY_LIMIT)
    cond
end


function accuracy_stop_check(status::State{<:xfgh_solution}, information, options)
    isnan(information.f_optimum) && (return false)
    
    vio = violationsSum(status.best_sol.g, status.best_sol.h, ε=options.h_tol)
    cond =  vio ≈ 0.0 && abs(status.best_sol.f - information.f_optimum) < options.f_tol
    cond && (status.termination_status_code = ACCURACY_LIMIT)
    cond
end

function accuracy_stop_check(status, information, options)
    return false
end

"""
    var_stop_check(status, information, options)

Check if the variance is close to zero in objective space.
"""
function var_stop_check(status, information, options)
    if typeof(status.best_sol.f) <: Array
        return false
    end
    cond =  var(map( s->s.f, status.population )) ≈ 0.0
    cond && (status.termination_status_code = OBJECTIVE_VARIANCE_LIMIT)
    cond
end

function diversity_stop_check(status, information, options)
    typeof(status.best_sol.f) <: Array && (return false)

    cond =  var(map( s->s.f, status.population )) ≈ 0.0
    cond && (status.termination_status_code = OBJECTIVE_VARIANCE_LIMIT)
    cond
end



function stop_check(status, information, options)
    accuracy_stop_check(status, information, options) ||
    var_stop_check(status, information, options)
end


"""
    diff_check(status, information, options; d = options.f_tol, p = 0.5)

Check the difference between best and worst objective function values in current
population (where at least %p of solution are feasible). Return `true` when such difference
is `<= d`, otherwise return `false`.


> Ref. Zielinski, K., & Laur, R. (n.d.). Stopping Criteria for Differential Evolution in
> Constrained Single-Objective Optimization. Studies in Computational Intelligence,
> 111–138. doi:10.1007/978-3-540-68830-3_4 (https://doi.org/10.1007/978-3-540-68830-3_4)
"""
function diff_check(status, information, options; d = options.f_tol, p = 0.3)
    !isnan(information.f_optimum) && (return false)

    mask = is_feasible.(status.population)
    p_feasibles = mean(mask)

    # not enough feasible samples?
    p_feasibles < p && (return false)

    fmin = minimum(s -> fval(s), status.population[mask])
    fmax = maximum(s -> fval(s), status.population[mask])

    cond = fmax - fmin <= d
    cond && (status.termination_status_code = OBJECTIVE_DIFFERENCE_LIMIT)
    cond
end



function default_stop_check(status, information, options)
    # budget termination
    call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)  ||
    time_stop_check(status, information, options) ||
    # accuracy
    accuracy_stop_check(status, information, options) ||
    # convergence
    stop_check(status, CheckConvergence(f_tol_abs = options.f_tol, x_tol=options.x_tol))
end
