struct Particle
    x::Vector
    f::Real
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

    return c / sum(fitness)
end

function getBest(Population::Array{Particle, 1}, searchType::Symbol)
    f_best = Population[1].f
    j = 1

    for i = 2:length(Population)
        if searchType == :minimize && f_best < Population[i].f
            f_best = Population[i].f
            j = i
        elseif searchType != :minimize && f_best > Population[i].f
            f_best = Population[i].f
            j = i
        end
    end

    return Population[j]
end


function eca(mfunc::Function,
                D::Int;
            η_max::Real= 4.0,
                K::Int = 3,
                N::Int = K * D,
        p_exploit::Real= 0.95,
            p_bin::Real= 0.3,
        max_evals::Int = 10000D,
      termination::Function = (x ->false),
      showResults::Bool = true,
       correctSol::Bool = true,
       searchType::Symbol=:minimize,
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

    for i in 1:N
        x = a + (b-a) .* rand(D)
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
    best = getBest(Population)
    
    convergence = []

    if saveConvergence != ""
        push!(convergence, best.f)
    end

    # start search
    while !stop

        # empty archive
        A   = Array{Particle, 1}([])
        
        I = randperm(N)

        p = nevals / max_evals

        # For each elements in Population
        for i in 1:N            
            x = Population[i].x

            # generate U masses
            if i <= N-K
                U_ids = I[i:K+i]
            else
                U_ids = I[1:K]
            end

            U = Population[U_ids]
            
            # generate center of mass
            c = center(U, searchType)

            # stepsize
            η = η_max * rand()

            # Ask if exploit process should be performed
            if p < p_exploit
                u = Population[rand(U_ids, 1)[1]].x

                # u-random-to-center/bin
                y = x + η * (c - u)

            else
                # center-to-best/bin
                y = x + η * (best.x - c)
            end

            for j = 1:D
                if rand() < p_bin
                    y[j] = x[j]
                     
                end
            end

            if correctSol
                y = correct(y, a, b)
            end

            f = func(y)
            sol = Particle(y, f)

            nevals += 1

            if Selection(Population[i], sol, searchType)
                push!(A, sol)

                # is it bestter?
                if Selection(best, sol, searchType)
                    best = sol
                end
            end
            
            stop = nevals >= max_evals
            if stop
                break
            end
        end

        t += 1

        replaceWorst!(Population, A, searchType)
        
        stop = stop || termination(Population)


        if saveConvergence != ""
            push!(convergence, best.f)
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
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        @printf("| best f.   = %e\n", best.f)
        println("=======================================")
    end

    return best.x, best.f
end