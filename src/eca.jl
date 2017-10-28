function center(neigbors, fitness)
    n, d = size(neigbors, 1, 2)
    c = zeros(Float64, d)

    fitness = 1.1maximum(fitness) - fitness


    for i = 1:n
        c += neigbors[i,:] * fitness[i]
    end

    return c / sum(fitness)
end

function replaceWorst!(population, fitness, A, f_A)
    f_wrost = sort(fitness, rev=true)

    l = 1
    for val in f_wrost[1:length(f_A)]
        j = find(x -> x == val, fitness)
        j = j[1]
        if fitness[j] < f_A[l]
            continue
        end
        population[j,:] = A[l]
        fitness[j] = f_A[l]
        l += 1
    end 
end

function eca(func::Function,
                D::Int,
            η_max::Real= 2.0,
                K::Int = 7,
                N::Int = K * D,
        max_evals::Int = 5000D,
      termination::Function = (x ->false),
      showResults::Bool = true,
           limits  = (-100., 100.))

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

            f_h = func(h)

            nevals += 1

            if f_h < fitness[i]
                push!(A,   h)
                push!(f_A, f_h)
            end
        end

        t += 1

        replaceWorst!(population, fitness, A, f_A)
        
        stop = nevals > max_evals || termination(fitness)
    end


    f_best = minimum(fitness)
    if showResults
        println("=======================================")
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        println("| best sol.   = ", f_best)
        println("| mean sol    = ", mean(fitness))
        println("| std. sol    = ", std(fitness))
        println("=======================================")
    end

    return population[find(x->x == f_best, fitness)[1], :], f_best
end