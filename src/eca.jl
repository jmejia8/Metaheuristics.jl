function replaceWorst!(Population::Array, A::Array, searchType::Symbol)
    n = length(A)
    if n == 0
        return
    end

    j = 1
    for i = randperm(length(Population))
        if Selection(Population[i], A[j], searchType)
            Population[i] = A[j]
            j += 1
        end
        if j>n
            return
        end
    end

end

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

function getMass(U::Array{xfg_indiv, 1}, searchType::Symbol)
    return getMass(U, searchType)
    
end

function getMass(U::Array{xfgh_indiv, 1}, searchType::Symbol)
    return getMass(U, searchType)
end

function center(U::Array, mass::Vector{Float64})
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

    func = fobj
    a, b = limits[1,:], limits[2,:]
    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end


    # population array
    Population = initializePop(func, N, D, a, b, initPopRand)

    # current evaluations
    nevals = N

    # stop condition
    stop = termination(Population)

    # current generation
    t = 0

    # best solution
    best = getBest(Population, searchType)

    convergence = []

    if saveConvergence != ""
        push!(convergence, [nevals best.f])
    end

    N_init = N
    p = nevals / max_evals

    if adaptive
        p_cr = rand(D)
    else
        p_cr = p_bin .* ones(D)
    end
    
    # start search
    while !stop
        I = randperm(N)

        adaptive && (Mcr_fail = zeros(D))
        
        # For each elements in Population
        for i in 1:N

            # current
            x = Population[i].x

            # generate U masses
            U = getU(Population, K, I, i, N)
            
            # generate center of mass
            c, u_worst, u_best = center(U, searchType)

            # stepsize
            η = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x

            # current-to-center/bin
            if p < p_exploit
                # u: worst element in U
                u = U[u_worst].x
                
                # current-to-center/bin
                y = x .+ η .* (c .- u)
            elseif p_exploit < 0
                y = x .+ (1-p^5)* η * (c .- u) .+ (p^5) * η * (best.x .- c)
            else
                # current-to-best/bin
                y = x .+ η .* (best.x .- c)
            end

            # binary crossover
            y, M_current = crossover(U[u_best].x, y, p_cr)

            y = correct(y, c, a, b)

            sol = generateChild(y, func(y))

            nevals += 1

            # replace worst element
            if Selection(Population[i], sol, searchType)
                Population[getWorstInd(Population, searchType)] = sol

                if Selection(best, sol, searchType)
                    best = sol
                end
            elseif adaptive
                Mcr_fail += M_current
            end
            
            stop = nevals >= max_evals
            if stop
                break
            end
        end

        if showIter
            @printf("| iter = %d \t nevals = %d \t f = %e\n", t, nevals, best.f)
            println("| ", best.x)
        end

        t += 1

        stop = stop || termination(Population)

        if stop
            break
        end

        adaptive && ( p_cr = adaptCrossover(p_cr, Mcr_fail/N) )


        if saveConvergence != ""
            push!(convergence, [nevals best.f])
        end

        p = nevals / max_evals
        
        if canResizePop
            # new size
            N = 2K .+ round(Int, (1 - p ) * (N_init .- 2K))

            if N < 2K
                N = 2K
            end

            Population = resizePop!(Population, N, K)
        end

    end

    if saveLast != ""
        o = []
        for i = 1:N
            push!(o, Population[i].x)
        end
        writecsv(saveLast, o)        
    end

    if saveConvergence != ""
        writecsv(saveConvergence, convergence)
    end


    if showResults
        println("===========[ ECA results ]=============")
        printResults(best, Population, t, nevals)
        println("=======================================")
    end

    if returnDetails
        return best, Population, t, nevals
    end

    return best.x, best.f
end
