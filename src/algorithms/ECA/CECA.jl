mutable struct CECA <: AbstractParameters
    η_max::Float64
    K::Int
    N::Int
    ε::Float64
end

"""
    CECA(;
        η_max = 2.0,
        K = 7,
        N = 0,
        ε = 0.0,
        information = Information(),
        options = Options()
    )

Parameters for the metaheuristic ECA: step-size `η_max`,`K` is number of vectors to
generate the center of mass, `N` is the population size.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], CECA())

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
function CECA(;
    η_max::Float64 = 2.0,
    K::Int = 7,
    N::Int = 0,
    ε::Float64 = 0.0,
    information = Information(),
    options = Options(),
)




    parameters = CECA(
        η_max,
        K,
        N,
        ε,
    )
    Algorithm(
        parameters,
        information = information,
        options = options,
    )

end


function center_ceca(U, max_fitness, ε)
    mass = getMass(U, max_fitness, ε)
    return center(U, mass),
    argworst(U), # worst
    argbest(U)  # best
end



function update_state!(
    status::State,
    parameters::CECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    K = parameters.K
    I = randperm(parameters.N)
    N = parameters.N
    population = status.population
    D = size(problem.bounds, 2)


    a = problem.bounds[1, :]
    b = problem.bounds[2, :]

    ε = 0.0
    feasible_solutions = findall( s->s.is_feasible, status.population )

    max_val = maximum( abs.(fvals(status)) )

    # For each elements in Population
    for i = 1:parameters.N
        p = status.f_calls / options.f_calls_limit

        if parameters.ε > 0
            ε = parameters.ε * (p - 1)^4
        end

        # current
        x = status.population[i].x

        # generate U masses
        U = getU(status.population, parameters.K, I, i, parameters.N, feasible_solutions)

        # generate center of mass
        c, u_worst, u_best = center_ceca(U, max_val, ε)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x


        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        y = x .+ η .* (c .- u)

        evo_boundary_repairer!(y, c, problem.bounds)

        sol = generateChild(y, problem.f(y))
        status.f_calls += 1

        # save new generated solution
        if is_better(sol, status.population[i])
            wi = argworst(status.population)
            status.population[wi] = sol

            if sol.is_feasible
                push!(feasible_solutions, wi)
            end
            


            if is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        else
            # stored but not used until replacement step
            push!(status.population, sol)
        end

        status.stop = stop_check(status, information, options)
        status.stop && break
    end
    
    if length(status.population) == N
        # no resize population
        return
    end

    # replacement step
    sort!(status.population, lt = is_better)
    deleteat!(status.population, N+1:length(status.population))
    


end


function initialize!(
    status::State,
    parameters::CECA,
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

    initialize!(problem, nothing, parameters, status, information, options)

end

function final_stage!(
    status::State,
    parameters::CECA,
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
        status.best_sol = get_pareto_front(status.population, is_better_eca)
    end
end

