function diffEvolution(func::Function, D::Int;
                        N::Int = 10D,
                        F::Real= 1.0,
                       CR::Real= 0.9,
                   CR_min::Real= CR,
                   CR_max::Real= CR,
                    F_min::Real=F,
                    F_max::Real=F,
                max_evals::Int = 10000D,
                 strategy::Symbol = :rand1,
              termination::Function = (x ->false),
               individual::DataType = xf_indiv,
              showResults::Bool = true,
                 saveLast::String = "",
          saveConvergence::String="",
                   limits  = [-100., 100.])

    if N < 5
       N = 5
       println("N increased to minimal value 5")
    end
    if CR < 0 || CR > 1
        CR = 0.5;
        println("CR should be from interval [0,1]; set to default value 0.5")
    end


    # bounds
    la, lb = limits[1,:], limits[2,:]
    if length(la) < D
        la = ones(D) * la[1]
        lb = ones(D) * lb[1]
    end

    # population array
    population = Array{individual, 1}([])

    # uniform initilization
    X = initializePop(N, D, la, lb)
    for i in 1:N
        x = X[i,:]
        push!(population, generateChild(individual, x, func(x)))
    end

    # current evalutations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0
    
    # best element ever
    best_ind = getBestInd(population)
    xBest = population[best_ind].x
    fBest = population[best_ind].f

    convergence = []

    if saveConvergence != ""
        push!(convergence, [nevals fBest])
    end

    # start search
    while !stop
        currentPop = copy(population)
        # For each elements in population

        # stepsize
        if F_min < F_max
            F = F_min + (F_max - F_min ) * rand()
        end

        if CR_min < CR_max
           CR = CR_min + (CR_max - CR_min) * rand()
        end

        for i in 1:N

            # select participats
            r1 = rand(1:N, 1)[1]
            while r1 == i
                r1 = rand(1:N, 1)[1]
            end

            r2 = rand(1:N, 1)[1]
            while r2 == i || r1 == r2
                r2 = rand(1:N, 1)[1]
            end

            r3 = rand(1:N, 1)[1]
            while r3 == i || r3 == r1 || r3 == r2
                r3 = rand(1:N, 1)[1]
            end

            x = currentPop[i].x
            a = currentPop[r1].x
            b = currentPop[r2].x
            c = currentPop[r3].x

            # strategy is selected here
            if strategy == :rand1
                # DE/rand/1
                u = a + F*(b - c)
            elseif strategy == :best1
                # DE/best/1
                u = xBest + F*(b - c)
            elseif strategy == :rand2
                # DE/rand/2

                r4 = rand(1:N, 1)[1]
                while r4 == i || r4 == r1 || r4 == r2 || r4 == r3
                    r4 = rand(1:N, 1)[1]
                end

                r5 = rand(1:N, 1)[1]
                while r5 == i || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4
                    r5 = rand(1:N, 1)[1]
                end

                d = currentPop[r4].x
                ee = currentPop[r5].x
                
                u = ee + F*(a - b + c - d)
            elseif strategy == :randToBest1
                # DE/rand-to-best/1
                u = x + F*(xBest - x + a - b)
            elseif strategy == :best2
                # DE/best/2
                r4 = rand(1:N, 1)[1]
                while r4 == i || r4 == r1 || r4 == r2 || r4 == r3 || r4 == best_ind
                    r4 = rand(1:N, 1)[1]
                end
                d = currentPop[r4].x
                u = xBest + F*(a - b + c - d)                
            else
                error("Unknown strategy $(strategy)")
            end

            # binomial crossover
            v = zeros(D)
            j_rand = rand(1:D,1)[1]

            # binomial crossover
            for j = 1:D
                if rand() < CR || j == j_rand
                    v[j] = u[j]

                    if v[j] < la[j]
                        v[j] = 2la[j] - v[j]
                    end
                    if v[j] > lb[j]
                        v[j] = 2lb[j] - v[j]
                    end
                else
                    v[j] = x[j]
                end
            end

            # instance child
            h = generateChild(individual, v, func(v))
            nevals += 1

            # select survivals
            if Selection(currentPop[i], h)
                population[i] = h

                if Selection(currentPop[best_ind], h)
                    best_ind = i
                end
            end

            stop = nevals >= max_evals
            if stop
                break
            end
        end

        t += 1

        xBest = population[best_ind].x
        fBest = population[best_ind].f

        if saveConvergence != ""
            push!(convergence, [nevals fBest])
        end
    
        # stop condition
        stop = stop || termination(population)
    end

    if saveLast != ""
        writecsv(saveLast, population)        
    end

    if saveConvergence != ""
        writecsv(saveConvergence, convergence)
    end

    if showResults
        println("============[ ED results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        println("| best sol.   = ", fBest)
        println("=======================================")
    end

    return xBest, fBest
end
