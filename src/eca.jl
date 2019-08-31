
function fitnessToMass(fitness::Vector{Float64}, searchType::Symbol)
    m = minimum(fitness)
    
    if m < 0
        fitness = 2abs(m) .+ fitness
    end

    if searchType == :minimize
        fitness = 2maximum(fitness) .- fitness
    end

    return fitness
end

function getMass(U::Array{xf_indiv, 1}, searchType::Symbol)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)
    
    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness, searchType)
end

function getMass(U::Array{xfg_indiv, 1}, searchType)
    return getMass(U, searchType)
    
end

function getMass(U::Array{xfgh_indiv, 1}, searchType)
    return getMass(U, searchType)
end

function center(U, mass)
    d = length(U[1].x)

    c = zeros(Float64, d)
    
    for i = 1:length(mass)
        c += mass[i] .* U[i].x
    end

    return c / sum(mass)
end

function center(U::Array, searchType::Symbol)
    n = length(U)

    mass = getMass(U, searchType)

    return center(U, mass), getWorstInd(U, searchType), getBestInd(U, searchType)
end

function getU(P::Array, K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N-K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j .+ 1]
    end

    return P[U_ids]
end

function crossover(x::Vector{Float64}, y::Vector{Float64}, p_cr::Vector{Float64})
    D = length(x)
    tmp2 = zeros(D)
    for j = 1:D
        if rand() < p_cr[j]
            y[j] = x[j]
            tmp2[j] += 1
        end
    end

    return y, tmp2
end

function adaptCrossover(p_cr::Vector{Float64}, M::Vector{Float64})
    p_best = p_cr[indmin(M)]

    for i = 1:length(p_cr)
        if M[i] > 0.3
            pn = abs(p_best .+ 0.3randn())
            if pn > 1.0
                pn = 1.0
            end

            if pn < 0.0
                pn = p_best
            end
            p_cr[i] = pn
        end
    end

    return p_cr
end

function resizePop!(P::Array, N_new::Int, K::Int)
    N = length(P)

    if N == N_new
        return P
    end

    f = zeros(N)
    for i = 1:N
        f[i] = P[i].f
    end

    ids = sortperm(P, lt=is_better)[1:N_new]
    return P[ids]
end

function update_state_eca!(problem, engine, parameters, status, information, options, iteration)
    K = parameters.K
    I = randperm(parameters.N)

    parameters.adaptive && (Mcr_fail = zeros(D))

    stop = false
    a = problem.bounds[1,:]
    b = problem.bounds[2,:]

    # For each elements in Population
    for i in 1:parameters.N
        p = status.f_calls / options.f_calls_limit

        # current
        x = status.population[i].x

        # generate U masses
        U = getU(status.population, parameters.K, I, i, parameters.N)
        
        # generate center of mass
        c, u_worst, u_best = center(U, :minimize)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        if p < parameters.p_exploit
            # u: worst element in U
            u = U[u_worst].x
            
            # current-to-center/bin
            y = x .+ η .* (c .- u)
        elseif parameters.p_exploit < 0
            y = x .+ (1-p^5)* η * (c .- u) .+ (p^5) * η * (status.best_sol.x .- c)
        else
            # current-to-best/bin
            y = x .+ η .* (status.best_sol.x .- c)
        end

        # binary crossover
        y, M_current = crossover(U[u_best].x, y, parameters.p_cr)

        y = correct(y, c, a, b)

        sol = generateChild(y, problem.f(y))

        status.f_calls += 1

        # replace worst element
        if engine.is_better(sol, status.population[i])
            status.population[getWorstInd(status.population, :minimize)] = sol
            if engine.is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        elseif parameters.adaptive
            Mcr_fail += M_current
        end
        
        status.stop = engine.stop_criteria(status, information, options) 
        status.stop && break
    end

    if status.stop
        return
    end

    if parameters.adaptive
        parameters.p_cr = adaptCrossover(parameters.p_cr, Mcr_fail/N)
    end

    # if saveConvergence != ""
    #     push!(convergence, [status.f_calls best.f])
    # end

    p = status.f_calls / options.f_calls_limit
    
    if parameters.resize_population
        # new size
        parameters.N = 2K .+ round(Int, (1 - p ) * (N_init .- 2K))

        if parameters.N < 2K
            parameters.N = 2K
        end

        status.population = resizePop!(status.population, parameters.N, K)
    end
end


function initialize_eca!(problem,engine,parameters,status,information,options)
    D = size(problem.bounds, 2)


    if parameters.N <= parameters.K
        parameters.N = parameters.K*D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&  @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    initialize!(problem,engine,parameters,status,information,options)

    N_init = parameters.N
    

    if parameters.adaptive
        parameters.p_cr = rand(D)
    else
        parameters.p_cr = parameters.p_bin .* ones(D)
    end
    
end

function final_stage_eca!(status, information, options)
    status.final_time = time()
    # if saveLast != ""
    #     o = []
    #     for i = 1:N
    #         push!(o, Population[i].x)
    #     end
    #     writecsv(saveLast, o)        
    # end

    # if saveConvergence != ""
    #     writecsv(saveConvergence, convergence)
    # end
end


