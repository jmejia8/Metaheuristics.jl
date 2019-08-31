function update_state_de!(problem, engine, parameters, status, information, options, iteration)
    population = status.population
    currentPop = copy(population)

    F = parameters.F
    CR = parameters.CR

    D = size(problem.bounds, 2)

    # stepsize
    if parameters.F_min < parameters.F_max
        F = parameters.F_min + (F_max - parameters.F_min ) * rand()
    end

    if parameters.CR_min < parameters.CR_max
       CR = parameters.CR_min + (parameters.CR_max - parameters.CR_min) * rand()
    end

    N = parameters.N
    strategy = parameters.strategy

    la = problem.bounds[1,:]
    lb = problem.bounds[2,:]

    D = length(la)

    xBest = status.best_sol.x

    for i in 1:N

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
            u = a + F*(b - c)
        elseif strategy == :best1
            # DE/best/1
            u = xBest + F*(b - c)
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
            
            u = ee + F*(a - b + c - d)
        elseif strategy == :randToBest1
            # DE/rand-to-best/1
            u = x + F*(xBest - x + a - b)
        elseif strategy == :best2
            # DE/best/2
            r4 = rand(1:N)
            while r4 == i || r4 == r1 || r4 == r2 || r4 == r3 || r4 == best_ind
                r4 = rand(1:N)
            end
            d = currentPop[r4].x
            u = xBest + F*(a - b + c - d)                
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

function initialize_de!(problem,engine,parameters,status,information,options)
    D = size(problem.bounds, 2)


    if parameters.N <= 5
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

end

function final_stage_de!(status, information, options)
    status.final_time = time()
end


function DE(f::Function, D::Int;
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
              showResults::Bool = true,
                 saveLast::String = "",
          saveConvergence::String="",
                   limits  = [-100., 100.])


    @warn "DE(f, D;...) function is deprecated. Please use: `optimize(f, bounds, DE())`"

    a, b = limits[1,:], limits[2,:]

    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')

    options = Options(f_calls_limit = max_evals, store_convergence = (saveConvergence != ""))
    method = DE()
    method.options.f_calls_limit = max_evals

    status = optimize(f, bounds, method)

    if true
        display(status)
    end


    return status.best_sol.x, status.best_sol.f
end

