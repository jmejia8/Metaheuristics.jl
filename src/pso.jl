function getBestPSO(fs::Array{Float64})
    fx = minimum(fs)

    return generateChild(zeros(2), fx)
end

function initialize_pso!(problem,engine,parameters,status,information,options)
    D = size(problem.bounds, 2)


    if parameters.N == 0
        parameters.N = 10*D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&  @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end



    initialize!(problem,engine,parameters,status,information,options)

    parameters.v = zeros(parameters.N, D)



    parameters.flock = status.population


end


@inline function velocity(x, v, pbest, gbest, parameters)
    r1 = parameters.C1*rand()
    r2 = parameters.C2*rand()
    parameters.ω*v + r2*(pbest - x) + r2*(gbest-x)
end

function update_state_pso!(problem, engine, parameters, status, information, options, iteration)
    xGBest = status.best_sol.x
    # For each elements in population
    for i in 1:parameters.N

        x = parameters.flock[i].x
        xPBest = status.population[i].x

        parameters.v[i,:] = velocity(x, parameters.v[i,:], xPBest, xGBest, parameters)
        x += parameters.v[i,:]

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

function pso(f::Function, D::Int;
                        N::Int = 10D,
                       C1::Real= 2.0,
                       C2::Real= 2.0,
                        ω::Real= 0.8,
                max_evals::Int = 10000D,
              termination::Function = (x ->false),
              showResults::Bool = true,
                   limits  = (-100., 100.))

    @warn "pso(f, D;...) function is deprecated. Please use: `optimize(f, bounds, PSO())`"

    a, b = limits

    if typeof(a) <: Real
        a = ones(D) * a
        b = ones(D) * b
    elseif length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')


    options = Options(f_calls_limit = max_evals)
    method = PSO(;N = N, C1 = C1, C2=C2, ω = ω)

    status = optimize(f, bounds, method)


    if true
        display(status)
    end

    return status.best_sol.x, status.best_sol.f
end
