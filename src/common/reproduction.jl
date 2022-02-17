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
        X[i,:] = ECA_operator(status.population,
                              parameters.K,
                              parameters.η_max;
                              i=i,
                             bounds = problem.bounds)
    end

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

