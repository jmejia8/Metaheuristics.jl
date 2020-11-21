function DE_mutation(population, i::Int, F::Float64, strategy::Symbol=:rand1, xBest=nothing)

    N = length(population)
    # select participats
    r1 = rand(1:N)
    while r1 == i
        r1 = rand(1:N)
    end

    r2 = rand(1:N)
    while r2 == i || r1 == r2
        r2 = rand(1:N)
    end

    r3 = rand(1:N)
    while r3 == i || r3 == r1 || r3 == r2
        r3 = rand(1:N)
    end

    x = population[i].x
    a = population[r1].x
    b = population[r2].x
    c = population[r3].x

    # strategy is selected here
    if strategy == :rand1
        # DE/rand/1
        u = a + F * (b - c)
    elseif strategy == :best1
        # DE/best/1
        u = xBest + F * (b - c)
    elseif strategy == :rand2
        # DE/rand/2

        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3
            r4 = rand(1:N)
        end

        r5 = rand(1:N)
        while r5 == i || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4
            r5 = rand(1:N)
        end

        d = population[r4].x
        ee = population[r5].x

        u = ee + F * (a - b + c - d)
    elseif strategy == :randToBest1
        # DE/rand-to-best/1
        u = x + F * (xBest - x + a - b)
    elseif strategy == :best2
        # DE/best/2
        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3 || r4 == best_ind
            r4 = rand(1:N)
        end
        d = population[r4].x
        u = xBest + F * (a - b + c - d)
    else
        @error("Unknown strategy $(strategy)")
    end

    return u
end

function DE_crossover(x, u, CR)
    D = length(x)
    # binomial crossover
    v = zeros(D)
    j_rand = rand(1:D)

    # binomial crossover
    for j = 1:D
        if rand() < CR || j == j_rand
            v[j] = u[j]
        else
            v[j] = x[j]
        end
    end

    return v
end

