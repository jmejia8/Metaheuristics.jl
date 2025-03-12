# structure with algorithm parameters
mutable struct MixedInteger{T} <: AbstractParameters
    base::T
    problem
end

include("utils.jl")

_num_to_vec(v::AbstractDict) = v
iscompatible(search_space::MixedSpace, algorithm::MixedInteger) = true

function MixedInteger(base; information = Information(), options = Options())
    parameters = MixedInteger(base.parameters, nothing)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end


function initialize!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    # introduce a helper objective function that translates vectors into dictionaries
    f(v) = begin
        d = vec_to_dict(v, problem.search_space)
        problem.f(d)
    end

    # create bounds from mixed
    bounds = _mixed_to_continous_space(problem.search_space)
    parameters.problem = Problem(f, bounds)


    status = initialize!(status, parameters.base, parameters.problem, information, options, args...; kargs...)

    #st = gen_initial_state(parameters.problem,parameters.base,information,options,status)

    problem.f_calls = parameters.problem.f_calls

    status
end


function update_state!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    update_state!(status, parameters.base, parameters.problem, information, options, args...; kargs...)
    problem.f_calls = parameters.problem.f_calls
end

function final_stage!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    
    # TODO: convert vectors in population into dictionaries
end


