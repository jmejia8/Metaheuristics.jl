mutable struct εDE <: AbstractDifferentialEvolution
    de::DE
    N::Int
    ε::Float64
    ε_0::Float64
    Tc::Int
    cp::Int
end

"""
    εDE(cp = 5, DE_kargs...)
    epsilonDE(cp = 5, DE_kargs...)

Parameters for ε Differential Evolution for constrained optimization.

See [DE](@ref) for more details about DE parameters (`DE_kargs`).

This implementation is not implementing the gradient-based repair method.
"""
function εDE(;cp = 5, kargs...)
    # get default Differential Evolution parameters
    de = DE(;kargs...)
    parameters = εDE(de.parameters, 0, 0.0, 0.0, 0, cp)
    info = de.information
    opts = de.options

    Algorithm(parameters, information=info, options=opts)
end

const epsilonDE = εDE

function is_better(a, b, parameters::εDE)
    ϕ_a = sum_violations(a)
    ϕ_b = sum_violations(b)
    ε = parameters.ε

    # ε level comparison
    if (ϕ_a <= ε && ϕ_b <= ε) || (ϕ_a == ϕ_b)
        return fval(a) < fval(b)
    end

    ϕ_a < ϕ_b 
end

ε_level_control_function(ε_0, t, Tc, cp) = t < Tc ? ε_0*(1 - t/Tc)^cp : 0.0

function initialize!(
        status,
        parameters::εDE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    status = initialize!(status, parameters.de, problem, information, options,args...; kargs...)
    elite_sols = sortperm(status.population, lt = is_better)

    θ = round(Int, 0.2length(elite_sols))
    s = rand(options.rng, status.population[elite_sols[1:θ]])
    
    parameters.ε_0 = sum_violations(s)
    parameters.Tc = round(Int, 0.2*options.iterations)
    parameters.N = parameters.de.N

    status


end

function update_state!(
        status,
        parameters::εDE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    # ε_level_control_function
    ε_0 = parameters.ε_0
    t = status.iteration
    Tc = parameters.Tc
    cp = parameters.cp
    parameters.ε = ε_level_control_function(ε_0, t, Tc, cp)

    # update_state!(status, parameters.de, args...; kargs...)
    new_vectors = reproduction(status, parameters.de, problem)

    # evaluate solutions
    new_solutions = create_solutions(new_vectors, problem,ε=options.h_tol)
    append!(status.population, new_solutions)

    # reduce population
    environmental_selection!(status.population, parameters)
    status.best_sol = get_best(status.population)


end
