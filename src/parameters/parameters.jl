function set_user_parameters!(algo::AbstractAlgorithm; kargs...)
    parameters = algo.parameters
    options = algo.options
    # information = algo.information
    T = typeof(parameters)
    for (field, value) in kargs
        # set algorithm parameters
        if ismutable(parameters) && hasfield(T, field)
            setproperty!(parameters, field, value)
        # set optional parameters
        elseif ismutable(options) && hasfield(Options, field)
            setproperty!(options, field, value)
        # elseif ismutable(information) && hasfield(Information, field)
        #     setproperty!(information, field, value)
        else
            @warn "Parameter `$field = $value` never used."
        end
    end
    algo
end

