###########################################
## DE
###########################################
function reproduction(status, parameters::DE, problem)
    @assert !isempty(status.population)

    N = parameters.N
    D = length(get_position(status.population[1]))

    strategy = parameters.strategy
    xBest = get_position(status.best_sol)
    population = status.population
    F = parameters.F
    CR = parameters.CR

    X = zeros(N,D)

    for i in 1:N
        x = get_position(population[i])
        u = DE_mutation(population, F, strategy, 1)
        v = DE_crossover(x, u, CR)
        evo_boundary_repairer!(v, xBest, problem.bounds)
        X[i,:] = v
    end

    X 
end

###########################################
## ECA
###########################################
function reproduction(status, parameters::ECA, problem)
    @assert !isempty(status.population)

    N = parameters.N
    D = length(get_position(status.population[1]))

    X = zeros(N,D)

    for i in 1:N
        X[i,:] = ECA_operator(population,
                              parameters.K,
                              parameters.η_max;
                              i=i,
                             bounds = problem.bounds)
    end

    X 
end

###########################################
## PSO
###########################################
function reproduction(status, parameters::PSO, problem)
    @assert !isempty(status.population)

    if isempty(parameters.flock)
        parameters.flock = status.population
    end

    N = parameters.N
    D = length(get_position(status.population[1]))

    if isempty(parameters.v)
        parameters.v = zeros(parameters.N, D)
    end

    X = zeros(N,D)
    xGBest = get_position(status.best_sol)

    for i in 1:N
        x = get_position(parameters.flock[i])
        xPBest = get_position(status.population[i])
        parameters.v[i, :] = velocity(x, parameters.v[i, :], xPBest, xGBest, parameters)
        x += parameters.v[i, :]
        reset_to_violated_bounds!(x, problem.bounds)
        X[i,:] = x
    end
    # set parameters.flock = status.population, out-side this function

    X 
end


###########################################
## CGSA reproduction
###########################################
function reproduction(status, parameters::CGSA, problem,options=nothing)
    @assert !isempty(status.population)
    if isnothing(options)
        @warn "Reproduction require options: reproduction(status,parameters,problem,options)"
        max_it = 500
    else
        max_it = options.iterations
    end
 
    # main parameters
    wMax = parameters.wMax
    N = parameters.N
    wMin = parameters.wMin
    iteration = status.iteration
    searchType = :minimize
    chaosIndex = parameters.chaosIndex	
    Rnorm = parameters.Rnorm
    Rpower = parameters.Rpower
    ElitistCheck = parameters.ElitistCheck

    # Velocity
    parameters.V = isempty(parameters.V) ? zeros(N,D) : parameters.V
    # Positions
    parameters.X = isempty(parameters.X) ? positions(status) : parameters.X
    # function values
    if isempty(parameters.fitness)
        parameters.fitness = fvals(status.population)
    end

    X = parameters.X
    V = parameters.V
    fitness = parameters.fitness

    P = status.population
    theBest = status.best_sol

    # iteration
    chValue = wMax-iteration*((wMax-wMin)/max_it)

    #Calculation of M. eq.14-20
    M = massCalculation(fitness,searchType)

    #Calculation of Gravitational constant. eq.13.
    G = Gconstant(iteration, max_it)

    if 1 <= chaosIndex <= 10
        G += chaos(chaosIndex,iteration,max_it,chValue)
    end


    #Calculation of accelaration in gravitational field. eq.7-10,21.
    a = Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)

    #Agent movement. eq.11-12
    X, V = move(X,a,V)

    for i = 1:N
        x = reset_to_violated_bounds!(X[i,:], problem.bounds)
        X[i,:] = x
    end

    parameters.X = X
    parameters.V = V

    X
end

###########################################
## generic GA reproduction
###########################################
function reproduction(status, parameters::AbstractNSGA, problem)
    @assert !isempty(status.population)

    I = randperm(parameters.N)
    J = randperm(parameters.N)
    Q = zeros(parameters.N, size(problem.bounds, 2))
    for i = 1:2:parameters.N

        pa = tournament_selection(status.population, I[i])
        pb = tournament_selection(status.population, J[i])

        c1, c2 = GA_reproduction(get_position(pa),
                                 get_position(pb),
                                 problem.bounds;
                                 η_cr = parameters.η_cr,
                                 p_cr = parameters.p_cr,
                                 η_m = parameters.η_m,
                                 p_m = parameters.p_m)
        Q[i,:] = c1
        Q[i+1,:] = c2
    end

    Q
end


###########################################
## NSGA3 reproduction
###########################################
function reproduction(status, parameters::NSGA3, problem)
    @assert !isempty(status.population)

    I = randperm(parameters.N)
    Q = zeros(parameters.N, size(problem.bounds, 2))
    for i = 1:parameters.N ÷ 2

        pa = status.population[I[2i-1]]
        pb = status.population[I[2i]]

        c1, c2 = GA_reproduction(get_position(pa),
                                 get_position(pb),
                                 problem.bounds;
                                 η_cr = parameters.η_cr,
                                 p_cr = parameters.p_cr,
                                 η_m = parameters.η_m,
                                 p_m = parameters.p_m)
        Q[i,:] = c1
        Q[i+1,:] = c2       
    end

    Q
end
