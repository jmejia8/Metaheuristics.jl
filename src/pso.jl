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
    stop = false

    # current generation
    t = 0
    
    i_min = indmin(fitness)
    xGBest= population[i_min,:]         
    fBest = fitness[i_min]

    velocity(x, v, pbest, gbest) = ω*v .+ C1*rand()*(pbest - x) .+ C2*rand()*(gbest-x)

    # start search
    while !stop
        # For each elements in population
        for i in 1:N

            x = population[i,:]
            xPBest = xPBests[i,:]

            v = velocity(x, vPop[i,:], xPBest, xGBest)
            x = x .+ v

            population[i,:] = x
            vPop[i,:] = v

            fx = func(x)
            nevals += 1

            if fx < fitness[i]
                xPBests[i,:]= x
                fitness[i]  = fx

                if fx < fBest
                    xGBest= x
                    fBest = fx
                end
            end

        end

        t += 1

        # stop condition
        stop = nevals > max_evals || termination(fitness)
    end


    if showResults
        println("============[ PSO results ]============")
        println("| Generations = $t")
        println("| Evals.      = ", nevals)
        println("| best sol.   = ", fBest)
        println("| mean sol    = ", mean(fitness))
        println("| std. sol    = ", std(fitness))
        println("=======================================")
    end

    return xGBest, fBest
end

# pso(x->sum(x.*x), 2;limits  = (-1., 1.))
