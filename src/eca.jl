struct Particle
    x::Vector{Float64}
    f::Float64
end

function Selection(fOld::Particle, fNew::Particle, searchType::Symbol)
    if searchType == :minimize
        return fNew.f < fOld.f
    end
    
    return fNew.f > fOld.f
end

function replaceWorst!(Population::Array{Particle, 1}, A::Array{Particle, 1}, searchType::Symbol)
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

function center(U::Array{Particle, 1}, searchType::Symbol)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)
    
    for i = 1:n
        fitness[i] = U[i].f
    end

    m = minimum(fitness)
    
    if m < 0
        fitness = 2abs(m) + fitness
    end

    if searchType == :minimize
        fitness = 2maximum(fitness) - fitness
    end

    c = zeros(Float64, d)
    for i = 1:n
        c += U[i].x * fitness[i]
    end

    return c / sum(fitness), indmin(fitness), indmax(fitness)
end

function getBest(Population::Array{Particle, 1}, searchType::Symbol)
    f_best = Population[1].f
    j = 1

    for i = 2:length(Population)
        if searchType == :minimize && f_best > Population[i].f
            f_best = Population[i].f
            j = i
        elseif searchType != :minimize && f_best < Population[i].f
            f_best = Population[i].f
            j = i
        end
    end

    return Population[j]
end

function getWorstInd(Population::Array{Particle, 1}, searchType::Symbol)
    f_worst = Population[1].f
    j = 1

    for i = 2:length(Population)
        if searchType == :minimize && f_worst < Population[i].f
            f_worst = Population[i].f
            j = i
        elseif searchType != :minimize && f_worst > Population[i].f
            f_worst = Population[i].f
            j = i
        end
    end

    return j
end

function getU(P::Array{Particle, 1}, K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N-K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j + 1]
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
            pn = abs(p_best + 0.3randn())
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

function resizePop!(P::Array{Particle, 1}, N_new::Int, K::Int)
    N = length(P)

    if N == N_new
        return P
    end

    f = zeros(N)
    for i = 1:N
        f[i] = P[i].f
    end

    ids = sortperm(f)[1:N_new]
    return P[ids]
end

function eca(mfunc::Function,
                D::Int;
            η_max::Real= 2,
                K::Int = 10,
                N::Int = 7*D,
        max_evals::Int = 10000D,
      termination::Function = (x ->false),
      showResults::Bool = true,
       correctSol::Bool = true,
       searchType::Symbol=:minimize,
       showIter::Bool = false,
         saveLast::String = "",
       saveConvergence::String="",
           limits  = [-100., 100.])

    func = mfunc
    a, b = limits[1,:], limits[2,:]
    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end


    Population = Array{Particle, 1}([])

    X = initializePop(N, D, a, b, :cheb)
    for i in 1:N
        x = X[i,:]
        f = func(x)
        push!(Population, Particle(x, f))
    end

    # current evaluations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solution
    best = getBest(Population, searchType)

    convergence = []

    if saveConvergence != ""
        push!(convergence, [nevals best.f])
    end

    N_init = N

    p_cr = rand(D)
    
    # start search
    while !stop
        I = randperm(N)

        p = nevals / max_evals

        Mcr_fail = zeros(D)
        
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
            y = x + η * (c - u)

            # binary crossover
            y, M_current = crossover(U[u_best].x, y, p_cr)

            y = correct(y, a, b, correctSol)

            f = func(y)
            sol = Particle(y, f)

            nevals += 1

            # replace worst element
            if Selection(Population[i], sol, searchType)
                Population[getWorstInd(Population, searchType)] = sol

                if Selection(best, sol, searchType)
                    best = sol
                end
            else
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

        p_cr = adaptCrossover(p_cr, Mcr_fail/N)


        if saveConvergence != ""
            push!(convergence, [nevals best.f])
        end

        # new size
        N = K + round(Int, (1-p)*(N_init - K))

        if N < K
            N = K
        end

        Population = resizePop!(Population, N, K)

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
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        @printf("| best f.     = %e\n", best.f)
        println("=======================================")
    end

    return best.x, best.f
end