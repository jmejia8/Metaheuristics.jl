function initialize!(
                status,
                parameters::AbstractParameters,
                problem,
                information,
                options,
                args...;
                kargs...
        )


    m = string(typeof(parameters))
    println("Method $m is calling the template function.")
    println("Redefine your own function using template")
    error("Method $m is calling the template function.")
    return State(0.0, zeros(0))
end


function update_state!(
        status,
        parameters::AbstractParameters,
        problem,
        information,
        options,
        args...;
        kargs...
)
        
    return
end



function final_stage(
        status,
        parameters::AbstractParameters,
        problem,
        information,
        options,
        args...;
        kargs...
)
    return
end

function stop_criteria!(
        status,
        parameters::AbstractParameters,
        problem,
        information,
        options,
        args...;
        kargs...
    )

        
    status.stop = call_limit_stop_check(status, information, options) ||
                  iteration_stop_check(status, information, options)  ||
                  time_stop_check(status, information, options)       ||
                  accuracy_stop_check(status, information, options)
    return
end

