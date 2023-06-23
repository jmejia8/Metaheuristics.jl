include("insert_random.jl")


"""
    DE_mutation(population, F = 1.0, strategy = :rand1)

Generate a `Vector` computed from population used in Differential Evolution.
Parameters: `F` is the stepsize, `strategy` can be one the following `:best1`, `:rand2`,
`:randToBest1` or `:best2`.
"""
function DE_mutation(population,
                    F::Float64 = 1.0,
                    strategy::Symbol=:rand1,
                    best_ind=0
    )

    N = length(population)

    i = rand(1:N)
    # select participats
    r1 = rand(1:N)
    while r1 == i
        r1 = rand(1:N)
    end

    r2 = rand(1:N)
    while r2 == i || r1 == r2
        r2 = rand(1:N)
    end


    x = get_position(population[i] )
    a = get_position(population[r1])
    b = get_position(population[r2])

    # strategy is selected here
    if strategy == :rand1
        # DE/rand/1
        return  x + F * (a - b)
    end

    r3 = rand(1:N)
    while r3 == i || r3 == r1 || r3 == r2
        r3 = rand(1:N)
    end

    c = get_position(population[r3])

    if strategy == :rand2
        # DE/rand/2

        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3
            r4 = rand(1:N)
        end

        r5 = rand(1:N)
        while r5 == i || r5 == r1 || r5 == r2 || r5 == r3 || r5 == r4
            r5 = rand(1:N)
        end

        d = get_position(population[r4])
        ee = get_position(population[r5])

        return ee + F * (a - b + c - d)
    end

    best_ind = best_ind == 0 ? argbest(population) : best_ind
    x_best = get_position(population[best_ind])
    if strategy == :best1
        # DE/best/1
        u = x_best + F * (b - c)
    elseif strategy == :randToBest1
        # DE/rand-to-best/1
        u = x + F * (x_best - x + a - b)
    elseif strategy == :best2
        # DE/best/2
        r4 = rand(1:N)
        while r4 == i || r4 == r1 || r4 == r2 || r4 == r3 || r4 == best_ind
            r4 = rand(1:N)
        end
        d = population[r4].x
        u = x_best + F * (a - b + c - d)
    else
        error("Unknown strategy " * string(strategy))
    end

    return u
end

