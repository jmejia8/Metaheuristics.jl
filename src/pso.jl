mutable struct PSO
    N::Int
    C1::Float64
    C2::Float64
    ω::Float64
    v::Array{Float64} # velocity
    flock::Array{xf_indiv}
end

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
        initialize! = initialize_pso!,
        update_state! = update_state_pso!,
        is_better = is_better,
        stop_criteria = stop_check,
        final_stage! = final_stage_pso!,
        information = information,
        options = options,
    )
end

function getBestPSO(fs::Array{Float64})
    fx = minimum(fs)

    return generateChild(zeros(2), fx)
end

function initialize_pso!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
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



    initialize!(problem, engine, parameters, status, information, options)

    parameters.v = zeros(parameters.N, D)



    parameters.flock = status.population


end


@inline function velocity(x, v, pbest, gbest, parameters)
    r1 = parameters.C1 * rand()
    r2 = parameters.C2 * rand()
    parameters.ω * v + r2 * (pbest - x) + r2 * (gbest - x)
end

function update_state_pso!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    iteration,
)
    xGBest = status.best_sol.x
    # For each elements in population
    for i = 1:parameters.N
        x = parameters.flock[i].x
        xPBest = status.population[i].x

        parameters.v[i, :] =
            velocity(x, parameters.v[i, :], xPBest, xGBest, parameters)
        x += parameters.v[i, :]

        sol = generateChild(x, problem.f(x))
        status.f_calls += 1

        if engine.is_better(sol, status.population[i])
            status.population[i] = sol

            if engine.is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        end

        parameters.flock[i] = sol

        # stop condition
        status.stop = engine.stop_criteria(status, information, options)
        status.stop && break
    end

end

function final_stage_pso!(status, information, options)
    status.final_time = time()

end
