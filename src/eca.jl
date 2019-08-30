
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
    println(status)
    Population = status.population
    K = parameters.K
    I = randperm(N)

    parameters.adaptive && (Mcr_fail = zeros(D))

    stop = false

    # For each elements in Population
    for i in 1:parameters.N

        # current
        x = Population[i].x

        # generate U masses
        U = getU(Population, parameters.K, I, i, parameters.N)
        
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
        if engine.is_better(Population[i], sol)
            Population[getWorstInd(Population, :minimize)] = sol

            if engine.is_better(status.best_sol, sol)
                status.best_sol = sol
            end
        elseif parameters.adaptive
            Mcr_fail += M_current
        end
        
        stop = stop_criteria(status, information, options) 
        stop && break
    end

    if showIter
        @printf("| iter = %d \t status.f_calls = %d \t f = %e\n", t, status.f_calls, best.f)
        println("| ", best.x)
    end

    t += 1

    stop = stop || stop_criteria(status, information, options) 

    if stop
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

        Population = resizePop!(Population, N, K)
    end
end

function initialize_eca!(problem,engine,parameters,status,information,options)
    a, b = problem.bounds[1,:], problem.bounds[2,:]
    D = length(a)

    # population array
    Population = initializePop(problem.f, parameters.N, D, a, b)

    # current evaluations
    status.f_calls = parameters.N

    # stop condition
    stop = engine.stop_criteria(status, information, options)

    # current generation
    t = 0

    # best solution
    best = getBest(Population, :minimize)
    status = State(best, Population)

    convergence = []

    if options.store_convergence
        status_tmp = deepcopy(status)
        empty!(status_tmp.convergence)

        push!(status.convergence, status_tmp)
    end

    N_init = parameters.N
    p = status.f_calls / options.f_calls_limit

    if parameters.adaptive
        parameters.p_cr = rand(D)
    else
        parameters.p_cr = parameters.p_bin .* ones(D)
    end
    
end

function final_stage_eca!(status, information, options)
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


function eca(fobj::Function,
                D::Int;
            η_max::Real  = 2,
                K::Int   = 7,
                N::Int   = K*D,
        p_exploit::Real  = 0.95,
            p_bin::Real  = 0.02,
        max_evals::Int   = 10000D,
      showResults::Bool  = true,
       correctSol::Bool  = true,
       searchType::Symbol=:minimize,
      initPopRand::Symbol=:uniform,
         showIter::Bool  = false,
         saveLast::String= "",
         adaptive::Bool  = false,
     canResizePop::Bool  = false,
      termination::Function   = (x ->false),
       saveConvergence::String="",
    returnDetails::Bool = false,
           limits  = [-100., 100.])

    @warn "This function is deprecated."

    f = fobj
    
    a, b = limits[1,:], limits[2,:]

    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')

    method = ECA()

    status = optimize(f, bounds, method)
    return status.best_sol.x, status.best_sol.f
end
