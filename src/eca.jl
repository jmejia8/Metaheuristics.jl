function center(neigbors::Matrix, fitness::Vector)
    n, d = size(neigbors, 1, 2)
    c = zeros(Float64, d)

    for i = 1:n
        c += neigbors[i,:] * fitness[i]
    end

    return c / sum(fitness)
end

function myFind(x, X)
    for i = 1:length(X)
        if x == X[i]
            return i
        end
    end

    return -1
end

function replaceWorst!(population::Matrix, fitness::Vector, A, f_A)
    f_wrost = sort(fitness)

    l = 1
    for val in f_wrost[1:length(f_A)]
        j = myFind(val, fitness)

        if fitness[j] > f_A[l]
            continue
        end

        population[j,:] = A[l]
        fitness[j] = f_A[l]
        l += 1
    end 
end

function correct(h, limits)
    a, b = limits
    return map(x-> begin 
                while(abs(x) > b)
                    x /= 2.0
                end
                return x
            end , h)
end

function eca(mfunc::Function,
                D::Int;
            η_max::Real= 2.0,
                K::Int = 7,
                N::Int = 2K * D,
        max_evals::Int = 10000D,
      termination::Function = (x ->false),
      showResults::Bool = true,
       correctSol::Bool = true,
         saveGens::Bool = false,
       searchType::Symbol=:minimize,
           limits  = (-100., 100.))

    func(x) = mfunc(x)
    if searchType == :minimize
        func(x) = 1.0 ./ (1 + mfunc(x))
    end

    a, b = limits

    population = a + (b - a) * rand(N, D)


    fitness = zeros(Float64, N)
    for i in 1:N            
        fitness[i] = func(population[i, :])
    end

    # current evalutations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solutions
    bestPerGen = []
    tmpBest = maximum(fitness)
    push!(bestPerGen, population[find(x->x == tmpBest, fitness)[1], :])

    # start search
    while !stop

        # empty archive
        A   = []
        f_A = Float64.([])
        
        # For each elements in population
        for i in 1:N            
            x = population[i, :]

            U_ids = rand(1:N, K)
            U = population[U_ids, :]

            η = η_max * rand()
            c = center(U, fitness[U_ids])
            u = population[rand(U_ids, 1)[1], :]

            h = x + η * (c - u)

            if correctSol
                h = correct(h, limits)
            end

            f_h = func(h)

            nevals += 1

            if f_h > fitness[i]
                push!(A,   h)
                push!(f_A, f_h)
            end
        end

        t += 1

        replaceWorst!(population, fitness, A, f_A)
        
        stop = nevals > max_evals || termination(fitness)

        if saveGens
            tmpBest = maximum(fitness)
            push!(bestPerGen, population[myFind(tmpBest, fitness), :])
        end
    end

    if saveGens
        writecsv("./solutions.csv", bestPerGen)        
    end

    fitness = -1.0 + 1.0 ./ fitness
    f_best = minimum(fitness)
    if showResults
        println("===========[ ECA results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        println("| best sol.   = ", f_best)
        println("| mean sol    = ", mean(fitness))
        println("| std. sol    = ", std(fitness))
        println("=======================================")
    end

    return population[myFind(f_best, fitness), :], f_best
end