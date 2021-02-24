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
`N` is the population size and `ω` controlls the inertia weight. 
# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO())
+=========== RESULT ==========+
| Iter.: 999
| f(x) = 3.23944e-48
| solution.x = [1.0698542573895642e-24, -1.4298101555926563e-24, -2.247029420442994e-25]
| f calls: 30000
| Total time: 0.4973 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO(N = 100, C1=1.5, C2=1.5, ω = 0.7))
+=========== RESULT ==========+
| Iter.: 299
| f(x) = 1.41505e-38
| solution.x = [2.161357427851024e-20, -1.1599444038307776e-19, 1.5122345732802047e-20]
| f calls: 30000
| Total time: 0.2128 s
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

function initialize!(
    status::State,
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



    initialize!(problem, nothing, parameters, status, information, options)

    parameters.v = zeros(parameters.N, D)



    parameters.flock = status.population


end


function update_state!(
    status::State,
    parameters::PSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
       )
    xGBest = status.best_sol.x
    # For each elements in population
    for i = 1:parameters.N
        x = parameters.flock[i].x
        xPBest = status.population[i].x

        parameters.v[i, :] =
        velocity(x, parameters.v[i, :], xPBest, xGBest, parameters)
        x = reset_to_violated_bounds!(x + parameters.v[i, :], problem.bounds)
        # x += parameters.v[i, :]

        sol = generateChild(x, problem.f(x), ε = options.h_tol)
        status.f_calls += 1

        if is_better(sol, status.population[i])
            status.population[i] = sol

            if is_better(sol, status.best_sol)
                status.best_sol = sol
                xGBest = status.best_sol.x
            end
        end

        parameters.flock[i] = sol

        # stop condition
        stop_criteria!(status, parameters, problem, information, options)
        status.stop && break
    end

end

function final_stage!(
    status::State,
    parameters::PSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()

end
