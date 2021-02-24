include("crossover_mutation.jl")

mutable struct DE <: AbstractParameters
    N::Int
    F::Float64
    CR::Float64
    CR_min::Float64
    CR_max::Float64
    F_min::Float64
    F_max::Float64
    strategy::Symbol
end

"""
    DE(;
        N  = 0,
        F  = 1.0,
        CR = 0.9,
        strategy = :rand1,
        information = Information(),
        options = Options()
    )

Parameters for Differential Evolution (DE) algorithm: step-size `F`,`CR` controlls the binomial
crossover, `N` is the population size. The parameter `trategy` is related to the variation
operator (`:rand1`, `:rand2`, `:best1`, `:best2`, `:randToBest1`).

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE())
+=========== RESULT ==========+
| Iter.: 437
| f(x) = 0
| solution.x = [0.0, 0.0, 0.0]
| f calls: 13134
| Total time: 0.3102 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE(N=50, F=1.5, CR=0.8))
+=========== RESULT ==========+
| Iter.: 599
| f(x) = 9.02214e-25
| solution.x = [-4.1003250484858545e-13, -6.090890160928905e-13, -6.025762626763004e-13]
| f calls: 30000
| Total time: 0.0616 s
+============================+
```

"""
function DE(;
    N::Int = 0,
    F = 1.0,
    CR = 0.9,
    CR_min = CR,
    CR_max = CR,
    F_min = F,
    F_max = F,
    strategy::Symbol = :rand1,
    information = Information(),
    options = Options(),
)


    parameters =
        DE(N, promote(F, CR, CR_min, CR_max, F_min, F_max)..., strategy)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )

end


function update_state!(
        status::State,
        parameters::DE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
    population = status.population
    currentPop = copy(population)

    F = parameters.F
    CR = parameters.CR


    # stepsize
    if parameters.F_min < parameters.F_max
        F = parameters.F_min + (F_max - parameters.F_min) * rand()
    end

    if parameters.CR_min < parameters.CR_max
        CR =
            parameters.CR_min + (parameters.CR_max - parameters.CR_min) * rand()
    end

    N = parameters.N
    strategy = parameters.strategy

    xBest = status.best_sol.x

    for i = 1:N

        if strategy == :randToBest1
            u = DE_mutation(currentPop, i, F, strategy, xBest)
        else
            u = DE_mutation(currentPop, i, F, strategy)
        end
        v = DE_crossover(currentPop[i].x, u, CR)

        # instance child
        v = evo_boundary_repairer!(v, xBest, problem.bounds)
        h = generateChild(v, problem.f(v),ε=options.h_tol)
        status.f_calls += 1

        # select survivals
        if is_better(h, currentPop[i])
            population[i] = h

            if is_better(h, status.best_sol)
                status.best_sol = h
                best_ind = i
            end
        end

        stop_criteria!(status, parameters, problem, information, options)
        if status.stop
            break
        end
    end




end

function initialize!(
        status::State,
        parameters::DE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
    D = size(problem.bounds, 2)



    if parameters.N <= 5
        parameters.N = 10 * D
    end

    if parameters.CR < 0 || parameters.CR > 1
        parameters.CR = 0.5
        options.debug &&
            @warn("CR should be from interval [0,1]; set to default value 0.5")
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
        parameters::DE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()
end
