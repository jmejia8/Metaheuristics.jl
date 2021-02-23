function initialize!(
        status::State,
        parameters::AbstractParameters,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
    return
end


function update_state!(
        status::State,
        parameters::AbstractParameters,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
        
    return
end



function final_stage(
        status::State,
        parameters::AbstractParameters,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
    return
end

function stop_criteria!(
        status::State,
        parameters::AbstractParameters,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

        
    status.stop = call_limit_stop_check(status, information, options) ||
                  iteration_stop_check(status, information, options)  ||
                  time_stop_check(status, information, options)
    return
end

