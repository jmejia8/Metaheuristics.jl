mutable struct Coevo <: AbstractParameters
    optimizers::Vector
    fitness::Vector
    migration_scheme
    problems::Vector{Problem}
    states::Vector{Status}
end


function initialize!(
    main_status,
    coevo::Coevo,
    main_problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...)

    coevo.problems = [Problem(fn, problem.bounds) for fn in coevo.fitness]

    for (status, optimizer, problem) in zip(coevo.states, coevo.optimizers, coevo.problems)
        initialize!(st, parameters, problem, information, options, args...; kargs...)
    end

end


function update_state!(
    main_status,
    parameters::Coevo,
    main_problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...)

    for (status, optimizer, problem) in zip(coevo.states, coevo.optimizers, coevo.problems)
        update_state!(status, parameters, problem, information, options, args...; kargs...)
    end

end


function final_stage!(
    main_status,
    parameters::Coevo,
    main_problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    for (status, optimizer, problem) in zip(coevo.states, coevo.optimizers, coevo.problems)
        final_stage!(status, parameters, problem, information, options, args...; kargs...)
    end
end
