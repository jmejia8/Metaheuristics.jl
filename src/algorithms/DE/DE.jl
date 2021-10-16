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
  iteration: 1000
    minimum: 0
  minimizer: [0.0, 0.0, 0.0]
    f calls: 30000
 total time: 0.0437 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE(N=50, F=1.5, CR=0.8))
+=========== RESULT ==========+
  iteration: 600
    minimum: 8.68798e-25
  minimizer: [3.2777877981303293e-13, 3.7650459509488005e-13, -7.871487597385812e-13]
    f calls: 30000
 total time: 0.0319 s
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
        status,
        parameters::DE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
)
    population = status.population
    current_pop = copy(population)

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

    xBest = get_position(status.best_sol)
    best_ind = 1

    for i in 1:N
        x = get_position(current_pop[i])

        u = DE_mutation(current_pop, F, strategy, best_ind)
        v = DE_crossover(x, u, CR)

        # instance child
        v = evo_boundary_repairer!(v, xBest, problem.bounds)
        h = create_solution(v, problem,ε=options.h_tol)

        # select survivals
        if is_better(h, current_pop[i])
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
        status,
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

    return gen_initial_state(problem,parameters,information,options,status)
end

function final_stage!(
        status,
        parameters::DE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()
end
