function call_limit_stop_check(status, information, options)
    cond =  status.f_calls >= options.f_calls_limit
    options.debug && cond && @info("Stopped since call_limit was met.")
    cond
end

function iteration_stop_check(status, information, options)
    cond =  status.iteration >= options.iterations
    options.debug && cond && @info("Stopped since iteration limit was met.")
    cond
end

function accuracy_stop_check(status, information, options)
    cond =  !isnan(information.f_optimum) && abs(status.best_sol.f - information.f_optimum) < options.f_tol
    options.debug && cond && @info("Stopped since accuracy was met.")
    cond
end

function var_stop_check(status, information, options)
    if typeof(status.best_sol.f) <: Array
        return false
    end
    cond =  var(map( s->s.f, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since varF was met.")
    cond
end

function diversity_stop_check(status, information, options)
    if typeof(status.best_sol.f) <: Array
        return false
    end

    cond =  var(map( s->s.f, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since varF was met.")
    cond
end



function stop_check(status, information, options)
    call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)    ||
    accuracy_stop_check(status, information, options) ||
    var_stop_check(status, information, options)
end
