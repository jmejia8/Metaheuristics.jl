const ε = 0.0001

struct Individual
    x::Vector
    f::Real
    g::Vector
    h::Vector
    νVal::Real
end


function G(g)
    i = g .<= 0.0
    
    gg = copy(g)
    gg[i] = 0.0
    return gg
end

function H(h)
    hh = abs.(h) - ε

    i = h .<= 0.0
    hh[i] = 0.0
    return hh
end

function ν(g, h)
    G_, H_ = G(g), H(h)
    m = length(G_) + length(H_)

    return (sum(G_) + sum(H_)) / m
end 

function is_feasible(ν_)
    return ν_ == 0.0
end

function Selection(xOld, xNew)
    new_feasible = is_feasible(xNew.νVal)
    old_feasible = is_feasible(xOld.νVal)



    if new_feasible && !old_feasible
        return true
    elseif new_feasible && old_feasible && xNew.f < xOld.f
        return true
    end


    if !new_feasible && !old_feasible && xNew.νVal < xOld.νVal
        return true
    end

    return  false
end

function replaceWorst!(Population::Array{Individual, 1}, A::Array{Individual, 1})
    n = length(A)
    if n == 0
        return
    end

    j = 1
    for i = randperm(length(Population))
        if Selection(Population[i], A[j])
            Population[i] = A[j]
            j += 1
        end
        if j>n
            return
        end
    end

end

function center(U::Array{Individual, 1}, searchType::Symbol)
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

function getBest(Population::Array{Individual, 1})
    ν_min = Population[1].νVal
    f_min = Population[1].f
    j = 1

    for i = 2:length(Population)
        if ν_min < Population[i].νVal 
            ν_min = Population[i].νVal
            f_min = Population[i].f
            j = i
        elseif ν_min == Population[i].νVal && f_min < Population[i].f
            ν_min = Population[i].νVal
            f_min = Population[i].f
            j = i
        end
    end

    return Population[j]
end

function ecaConstrained(
            mfunc::Function,
                D::Int,
                D_g::Int,
                D_h::Int;
            η_max::Real= 2.0,
                K::Int = 7,
                N::Int = 2K * D,
        max_evals::Int = 10000D,
      termination::Function = (x ->false),
      showResults::Bool = true,
       correctSol::Bool = true,
         saveGens::Bool = false,
     negativeVals::Bool = false,
       searchType::Symbol=:minimize,
       saveConvergence::Bool=false,
           limits  = [-100., 100.])

    func = mfunc

    a, b = limits[1,:], limits[2,:]
    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end


    Population = Array{Individual, 1}([])

    for i in 1:N
        x     = a + (b-a) .* rand(D)
        f,g,h = func(x)
        νVal  = ν(g,h)
        
        push!(Population, Individual(x, f, g, h, νVal))
    end

    # current evaluations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solutions
    # bestPerGen = []
    # convergence = []
    # tmpBest = maximum(f)

    if saveGens
        push!(bestPerGen, Population[find(x->x == tmpBest, f)[1], :])
    end

    if saveConvergence
        push!(convergence, tmpBest)
    end

    # start search
    while !stop

        # empty archive
        A   = Array{Individual, 1}([])
        
        # For each elements in Population
        for i in 1:N            
            x = Population[i].x

            U_ids = randperm(N)[1:K]
            U = Population[U_ids]

            η = η_max * rand()
            c = center(U, searchType)
            u = Population[rand(U_ids, 1)[1]].x

            y = x + η * (c - u)

            if correctSol
                y = correct(y, limits)
            end

            f, g, h = func(y)
            sol = Individual(y, f, g, h, abs(ν(g, h)))

            nevals += 1

            if Selection(Population[i], sol)
                push!(A, sol)
            end
        end

        t += 1

        replaceWorst!(Population, A)
        
        stop = nevals > max_evals || termination(Population)

        # if saveGens
        #     tmpBest = maximum(f)
        #     push!(bestPerGen, Population[myFind(tmpBest, f), :])
        # end

        # if saveConvergence
        #     tmpBest = maximum(f)
        #     push!(convergence, tmpBest)
        # end
    end

    # if saveGens
    #     writecsv("./solutions.csv", bestPerGen)        
    # end

    # if saveConvergence
    #     writecsv("./convergence.csv", convergence)
    # end

    best = getBest(Population)

    if showResults
        println("===========[ ECA results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", nevals)
        println("| best f.   = ", best.f)
        println("| best ν.   = ", best.νVal)
        println("=======================================")
    end

    return best
end