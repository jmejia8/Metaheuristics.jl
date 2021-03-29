include("mass.jl")
include("adaptive_parameters.jl")
include("center_of_mass.jl")

mutable struct ECA <: AbstractParameters
    η_max::Float64
    K::Int
    N::Int
    N_init::Int
    p_exploit::Float64
    p_bin::Float64
    p_cr::Array{Float64}
    adaptive::Bool
    resize_population::Bool
end

"""
    ECA(;
        η_max = 2.0,
        K = 7,
        N = 0,
        N_init = N,
        p_exploit = 0.95,
        p_bin = 0.02,
        p_cr = Float64[],
        adaptive = false,
        resize_population = false,
        information = Information(),
        options = Options()
    )

Parameters for the metaheuristic ECA: step-size `η_max`,`K` is number of vectors to
generate the center of mass, `N` is the population size.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], ECA())

+=========== RESULT ==========+
| Iter.: 1021
| f(x) = 1.68681e-163
| solution.x = [2.5517634463667404e-82, -2.9182760041942484e-82, -1.3565584801935802e-82]
| f calls: 21454
| Total time: 0.0894 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], ECA(N = 10, η_max = 1.0, K = 3))
+=========== RESULT ==========+
| Iter.: 1506
| f(x) = 0.000172391
| solution.x = [-6.340714627875324e-5, -0.004127226953894587, 0.012464071313908906]
| f calls: 15069
| Total time: 0.0531 s
+============================+
```

"""
function ECA(;
    η_max::Float64 = 2.0,
    K::Int = 7,
    N::Int = 0,
    N_init::Int = N,
    p_exploit::Float64 = 0.95,
    p_bin::Float64 = 0.02,
    p_cr::Array{Float64} = Float64[],
    adaptive::Bool = false,
    resize_population::Bool = false,
    information = Information(),
    options = Options(),
)


    N_init = N


    parameters = ECA(
        η_max,
        K,
        N,
        N_init,
        p_exploit,
        p_bin,
        p_cr,
        adaptive,
        resize_population,
    )
    Algorithm(
        parameters,
        information = information,
        options = options,
    )

end


function update_state!(
    status,
    parameters::ECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    K = parameters.K
    I = randperm(parameters.N)
    D = size(problem.bounds, 2)
    is_multiobjective = typeof(status.population[1]) == xFgh_indiv

    parameters.adaptive && (Mcr_fail = zeros(D))

    stop = false
    a = problem.bounds[1, :]
    b = problem.bounds[2, :]

    if is_multiobjective
        shuffle!(status.population)
    end

    # For each elements in Population
    for i = 1:parameters.N
        p = status.f_calls / options.f_calls_limit


        # current
        x = status.population[i].x

        # generate U masses
        U = getU(status.population, parameters.K, I, i, parameters.N)

        # generate center of mass
        c, u_worst, u_best = center(U)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        if p < parameters.p_exploit
            # u: worst element in U
            u = U[u_worst].x

            # current-to-center/bin
            y = x .+ η .* (c .- u)
        elseif parameters.p_exploit < 0
            y =
                x .+ (1 - p^5) * η * (c .- u) .+
                (p^5) * η * (status.best_sol.x .- c)
        else
            # current-to-best/bin
            y = x .+ η .* (status.best_sol.x .- c)
        end

        # binary crossover
        y, M_current = crossover(U[u_best].x, y, parameters.p_cr)

        evo_boundary_repairer!(y, c, problem.bounds)

        sol = generateChild(y, problem.f(y))

        status.f_calls += 1

        # replace worst element
        if is_better(sol, status.population[i])
            wi = argworst(status.population)
            status.population[wi] = sol
            if is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        else
            if is_multiobjective && !is_better(status.population[i], sol)
                push!(status.population, sol)
            end
            parameters.adaptive && (Mcr_fail += M_current)
        end

        stop_criteria!(status, parameters, problem, information, options)
        status.stop && break
    end

    if status.stop
        return
    end

    if parameters.adaptive
        parameters.p_cr =
            adaptCrossover(parameters.p_cr, Mcr_fail / parameters.N)
    end

    if is_multiobjective && length(status.population) > parameters.N
        status.population = truncate_population!(status.population, parameters.N, (w,z) -> is_better(w, z))
    end

    p = status.f_calls / options.f_calls_limit

    if parameters.resize_population
        K = parameters.K

        # new size
        parameters.N = 2K .+ round(Int, (1 - p) * (parameters.N_init .- 2K))

        if parameters.N < 2K
            parameters.N = 2K
        end

        status.population = resizePop!(
            status.population,
            parameters.N,
            K,
            (a, b) -> is_better(a, b),
        )
    end
end


function initialize!(
    parameters::ECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    D = size(problem.bounds, 2)


    if parameters.N <= parameters.K
        parameters.N = parameters.K * D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
            @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    # initialize!(problem, nothing, parameters, status, information, options)
    status = gen_initial_state(problem,parameters,information,options)

    N_init = parameters.N


    if parameters.adaptive
        parameters.p_cr = rand(D)
    else
        parameters.p_cr = parameters.p_bin .* ones(D)
    end

    status

end

function final_stage!(
    status,
    parameters::ECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()

    # compute Pareto front if it is a multiobjective problem
    if typeof(status.population[1].f) <: Array
        options.debug && @info "Computing Pareto front..."
        status.best_sol = get_pareto_front(status.population, is_better)
    end
end


function is_better_eca(
    New::xFgh_indiv,
    Old::xFgh_indiv;
    searchType = :minimize,
    ε = 0.0,
)

    old_vio = Old.sum_violations
    new_vio = New.sum_violations

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    if searchType == :minimize
        for i in 1:length(Old.f)
            if Old.f[i] < New.f[i]
                return false
            end
        end
        return true
    end

    for i in 1:length(Old.f)
        if Old.f[i] > New.f[i]
            return false
        end
    end

    return true
end
