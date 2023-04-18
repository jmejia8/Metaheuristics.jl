#=
function stop_criteria!(
        status,
        parameters::Union{ECA, DE, PSO, CGSA},
        problem,
        information,
        options,
        args...;
        kargs...
    )
    status.stop && (return)

    status.stop = diff_check(status, information, options)
end
=#
