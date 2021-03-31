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

function time_stop_check(status, information, options)
    cond =  status.overall_time >= options.time_limit
    options.debug && cond && @info("Stopped since time limit was met.")
    cond
end

function accuracy_stop_check(status::State{xf_indiv}, information, options)
    cond =  !isnan(information.f_optimum) && abs(status.best_sol.f - information.f_optimum) < options.f_tol
    options.debug && cond && @info("Stopped since accuracy was met.")
    cond
end


function accuracy_stop_check(status::State{xfgh_indiv}, information, options)
    if isnan(information.f_optimum)
        return false
    end
    
    vio = violationsSum(status.best_sol.g, status.best_sol.h, ε=options.h_tol)
    cond =  vio ≈ 0.0 && abs(status.best_sol.f - information.f_optimum) < options.f_tol

    options.debug && cond && @info("Stopped since accuracy and feasibility was met.")
    cond
end

function accuracy_stop_check(status::State{xFgh_indiv}, information, options)
    return false
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
    accuracy_stop_check(status, information, options) ||
    var_stop_check(status, information, options)
end
