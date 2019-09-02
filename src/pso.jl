function getBestPSO(fs::Array{Float64})
    fx = minimum(fs)

    return generateChild(zeros(2), fx)
end


function velocity(x, v, pbest, gbest, parameters)
    parameters.ω*v .+ parameters.C1*rand()*(pbest - x) + parameters.C2*rand()*(gbest-x)
end

function update_state_pso!(problem, engine, parameters, status, information, options, iteration)
    D = size(problem.bounds, 2)


    if parameters.N <= parameters.K
        parameters.N = parameters.K*D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&  @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

end

function update_state_pso!(problem, engine, parameters, status, information, options, iteration)
    xGBest = status.best_sol.x
    # For each elements in population
    for i in 1:parameters.N

        x = status.population[i]
        xPBest = parameters.pbest[i,:]

        parameters.v[i] = velocity(x, vPop[i,:], xPBest, xGBest, parameters)
        x = x + v

        sol = generateChild(x, problem.f(x))
        status.f_calls += 1

        if engine.is_better(sol, status.population[i])
            parameters.pbest[i,:] = x

            if engine.is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        end

        status.population[i] = sol
        
        # stop condition
        status.stop = engine.stop_criteria(status, information, options) 
        status.stop && break
    end

end

function pso(func::Function, D::Int;
                        N::Int = 10D,
                       C1::Real= 2.0,
                       C2::Real= 2.0,
                        ω::Real= 0.8,
                max_evals::Int = 10000D,
              termination::Function = (x ->false),
              showResults::Bool = true,
                   limits  = (-100., 100.))

    a, b = limits

    population = a .+ (b .- a) * rand(N, D)
    vPop = randn(N, D)
   
    fitness = zeros(Real, N)
    for i in 1:N            
        fitness[i] = func(population[i, :])
    end

    xPBests = copy(population)

    # current evaluationsVu
    nevals = N

    # stop condition
    stop = termination(fitness)

    # current generation
    t = 0
    
    i_min = indmin(fitness)
    xGBest= population[i_min,:]         
    fBest = fitness[i_min]

    velocity(x, v, pbest, gbest) = ω*v .+ C1*rand()*(pbest - x) .+ C2*rand()*(gbest-x)

    # start search
    while !stop

    end


    if showResults
        display(status)
    end

    return xGBest, fBest
end

# pso(x->sum(x.*x), 2;limits  = (-1., 1.))
