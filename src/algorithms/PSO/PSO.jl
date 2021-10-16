include("velocity.jl")


mutable struct PSO <: AbstractParameters
    N::Int
    C1::Float64
    C2::Float64
    ω::Float64
    v::Array{Float64} # velocity
    flock::Array
end


"""
    PSO(;
        N  = 0,
        C1 = 2.0,
        C2 = 2.0,
        ω  = 0.8,
        information = Information(),
        options = Options()
    )

Parameters for Particle Swarm Optimization (PSO) algorithm: learning rates `C1` and `C2`,
`N` is the population size and `ω` controls the inertia weight. 
# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO())
+=========== RESULT ==========+
  iteration: 1000
    minimum: 1.40522e-49
  minimizer: [3.0325415595139883e-25, 1.9862212295897505e-25, 9.543772256546461e-26]
    f calls: 30000
 total time: 0.1558 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO(N = 100, C1=1.5, C2=1.5, ω = 0.7))
+=========== RESULT ==========+
  iteration: 300
    minimum: 2.46164e-39
  minimizer: [-3.055334698085433e-20, -8.666986835846171e-21, -3.8118413472544027e-20]
    f calls: 30000
 total time: 0.1365 s
+============================+
```

"""
function PSO(;
    N::Int = 0,
    C1 = 2.0,
    C2 = 2.0,
    ω = 0.8,
    v = Float64[],
    flock = xf_indiv[],
    information = Information(),
    options = Options(),
)

parameters = PSO(N, promote(Float64(C1), C2, ω)..., v, flock)

Algorithm(
    parameters,
    information = information,
    options = options,
)
end

function update_state!(
        status,
        parameters::PSO,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    xGBest = get_position(status.best_sol)

    # For each elements in population
    for i in 1:parameters.N
        x = get_position(parameters.flock[i])
        xPBest = get_position(status.population[i])

        parameters.v[i, :] = velocity(x, parameters.v[i, :], xPBest, xGBest, parameters)
        x += parameters.v[i, :]
        reset_to_violated_bounds!(x, problem.bounds)

        sol = create_solution(x, problem, ε = options.h_tol)

        if is_better(sol, status.population[i])
            status.population[i] = sol

            if is_better(sol, status.best_sol)
                status.best_sol = sol
                xGBest = get_position(status.best_sol)
            end
        end

        parameters.flock[i] = sol

        # stop condition
        stop_criteria!(status, parameters, problem, information, options)
        status.stop && break
    end

end

function initialize!(
    status,
    parameters::PSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
       )
    D = size(problem.bounds, 2)


    if parameters.N == 0
        parameters.N = 10 * D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
        @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end



    status = gen_initial_state(problem,parameters,information,options,status)

    parameters.v = zeros(parameters.N, D)



    parameters.flock = status.population

    status

end

function final_stage!(
    status,
    parameters::PSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()

end
