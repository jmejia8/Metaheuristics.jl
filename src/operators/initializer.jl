abstract type AbstractInitializer end

function initialize!(
        status,
        parameters::AbstractInitializer,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    gen_initial_state(problem,parameters,information,options, status)
end

