function diffEvolution(func::Function, D::Int;
                        N::Int = 10D,
                        F::Real= 2.0,
                       CR::Real= 0.9,
                max_evals::Int = 10000D,
                 strategy::Symbol = :rand1,
              termination::Function = (x ->false),
              showResults::Bool = true,
                   limits  = (-100., 100.))

    if N < 5
       N = 5
       println("N increased to minimal value 5")
    end
    if CR < 0 || CR > 1
        CR = 0.5;
        println("CR should be from interval [0,1]; set to default value 0.5")
    end

    a, b = limits

    population = a + (b - a) * rand(N, D)

    fitness = zeros(Real, N)
    for i in 1:N            
        fitness[i] = func(population[i, :])
    end

    # current evalutations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0
    
    i_min = indmin(fitness)
    xBest = population[i_min,:]
    fBest = fitness[i_min]

    # start search
    while !stop
        x = copy(population)
        # For each elements in population
        for i in 1:N
            k = randperm(N)

            if strategy == :rand1
                # DE/rand/1
                u = x[k[3],:] + F*(x[k[1],:] - x[k[2],:])
            elseif strategy == :best1
                # DE/best/1
                u = xBest + F*(x[k[1],:] - x[k[2],:])
            elseif strategy == :rand2
                # DE/rand/2
                u = x[k[5],:] + F*(x[k[1],:] - x[k[2],:] + x[k[3],:] - x[k[4],:])
            elseif strategy == :randToBest1
                # DE/rand-to-best/1
                u = x[i,:] + F*(xBest - x[i,:] + x[k[1],:] - x[k[2],:])
            elseif strategy == :best2
                # DE/best/2
                u = xBest + F*(x[k[1],:] - x[k[2],:] + x[k[3],:] - x[k[4],:])                
            else
                error("Unknown strategy $(strategy)")
            end

            
            # binomial crossover
            v = x[i,:]
            r = rand(1:D,1)[1]

            I = rand(D) .< CR
            v[I] = u[I]
            v[r] = u[r]

            # Correct sol.
            I = .!(a .< v .< b)
            v[ I ] = a + (b-a) * rand(sum(I))

            fv = func(v)
            nevals += 1
            if fv < fitness[i]
                population[i,:] = v
                fitness[i] = fv
            end

            stop = nevals >= max_evals
            if stop
                break
            end
        end

        t += 1

        i_min = indmin(fitness)
        xBest = population[i_min,:]
        fBest = fitness[i_min]
    
        # stop condition
        stop = stop || termination(fitness)
    end


    if showResults
        println("============[ ED results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        println("| best sol.   = ", fBest)
        println("| mean sol    = ", mean(fitness))
        println("| std. sol    = ", std(fitness))
        println("=======================================")
    end

    return xBest, fBest
end