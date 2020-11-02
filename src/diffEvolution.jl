mutable struct DE
    N::Int
    F::Float64
    CR::Float64
    CR_min::Float64
    CR_max::Float64
    F_min::Float64
    F_max::Float64
    strategy::Symbol
end

"""
    DE(;
        N  = 0,
        F  = 1.0,
        CR = 0.9,
        strategy = :rand1,
        information = Information(),
        options = Options()
    )

Parameters for Differential Evolution (DE) algorithm: step-size `F`,`CR` controlls the binomial
crossover, `N` is the population size. The parameter `trategy` is related to the variation
operator (`:rand1`, `:rand2`, `:best1`, `:best2`, `:randToBest1`).

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE())
+=========== RESULT ==========+
| Iter.: 437
| f(x) = 0
| solution.x = [0.0, 0.0, 0.0]
| f calls: 13134
| Total time: 0.3102 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE(N=50, F=1.5, CR=0.8))
+=========== RESULT ==========+
| Iter.: 599
| f(x) = 9.02214e-25
| solution.x = [-4.1003250484858545e-13, -6.090890160928905e-13, -6.025762626763004e-13]
| f calls: 30000
| Total time: 0.0616 s
+============================+
```

"""
function DE(;
    N::Int = 0,
    F = 1.0,
    CR = 0.9,
    CR_min = CR,
    CR_max = CR,
    F_min = F,
    F_max = F,
    strategy::Symbol = :rand1,
    information = Information(),
    options = Options(),
)


    parameters =
        DE(N, promote(F, CR, CR_min, CR_max, F_min, F_max)..., strategy)

    Algorithm(
        parameters,
        initialize! = initialize_de!,
        update_state! = update_state_de!,
        is_better = is_better,
        stop_criteria = stop_check,
        final_stage! = final_stage_de!,
        information = information,
        options = options,
    )

end


function update_state_de!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    iteration,
)
    population = status.population
    currentPop = copy(population)

    F = parameters.F
    CR = parameters.CR

    D = size(problem.bounds, 2)

    # stepsize
    if parameters.F_min < parameters.F_max
        F = parameters.F_min + (F_max - parameters.F_min) * rand()
    end

    if parameters.CR_min < parameters.CR_max
        CR =
            parameters.CR_min + (parameters.CR_max - parameters.CR_min) * rand()
    end

    N = parameters.N
    strategy = parameters.strategy

    la = problem.bounds[1, :]
    lb = problem.bounds[2, :]

    D = length(la)

    xBest = status.best_sol.x

    for i = 1:N

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

        x = currentPop[i].x
        a = currentPop[r1].x
        b = currentPop[r2].x
        c = currentPop[r3].x

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

            d = currentPop[r4].x
            ee = currentPop[r5].x

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
            d = currentPop[r4].x
            u = xBest + F * (a - b + c - d)
        else
            @error("Unknown strategy $(strategy)")
        end

        # binomial crossover
        v = zeros(D)
        j_rand = rand(1:D)

        # binomial crossover
        for j = 1:D
            if rand() < CR || j == j_rand
                v[j] = u[j]

                if v[j] < la[j]
                    v[j] = la[j]
                elseif v[j] > lb[j]
                    v[j] = lb[j]
                end
            else
                v[j] = x[j]
            end
        end

        # instance child
        h = generateChild(v, problem.f(v))
        status.f_calls += 1

        # select survivals
        if engine.is_better(h, currentPop[i])
            population[i] = h

            if engine.is_better(h, status.best_sol)
                status.best_sol = h
                best_ind = i
            end
        end

        status.stop = engine.stop_criteria(status, information, options)
        if status.stop
            break
        end
    end




end

function initialize_de!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
)
    D = size(problem.bounds, 2)



    if parameters.N <= 5
        parameters.N = 10 * D
    end

    if parameters.CR < 0 || parameters.CR > 1
        parameters.CR = 0.5
        options.debug &&
            @warn("CR should be from interval [0,1]; set to default value 0.5")
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

end

function final_stage_de!(status, information, options)
    status.final_time = time()
end
