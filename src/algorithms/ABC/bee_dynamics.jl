mutable struct Bee{T}
    sol::T
    fit::Float64
    t::Int
end


Bees = Array{Bee, 1}

Bee(solution) = Bee(solution, fit(solution.f), 0)


"""
    get_position(bee)

Get the position vector of a bee when optimize using ABC algorithm.
"""
get_position(bee::Bee) = bee.sol.x

"""
    fval(solution)

Get the fitness of a bee when optimize using ABC algorithm.
"""
fval(bee::Bee) = bee.sol.f

minimum(st::State{Bee{xf_indiv}}) = st.best_sol.sol.f

@inline function updateFit!(bee::Bee)
    bee.fit = fit(bee.sol.f)
end

@inline function updateFit!(bees::Bees)
    for bee in bees
        updateFit!(bee)
    end
end

function fit(fx)
    if fx >= 0.0
        return 1.0/(1.0+fx)
    end
    
    1.0 + abs(fx)
end

function updateBee!(bee, bee2, f, bounds)
    D = length(bee.sol.x)
    ϕ = -1.0 + 2.0rand()

    v = ϕ*(bee.sol.x - bee2.sol.x)

    x_new = bee.sol.x + v

    fx = f(x_new)

    if fx < bee.sol.f
        bee.sol.x = x_new
        bee.sol.f = fx
        bee.fit = fit(fx)
        bee.t = 0
    else
        bee.t += 1
    end
end

function getK(i, N)
    return rand(union(1:i-1, i+1:N))
end

function employedPhase!(bees, f, Ne, bounds)
    N = length(bees)
    for i in randperm(N)[1:Ne]
        updateBee!(bees[i], bees[getK(i, N)], f, bounds)
    end
end


function roulettSelect(bees, sum_f)
    r = rand()
    rs = 0.0
    for i in 1:length(bees)
        rs += bees[i].fit / sum_f
        if r <= rs
            return i
        end
    end

    return length(bees)
end

function outlookerPhase!(bees, f, No::Int, bounds)
    N = length(bees)
    sum_f = sum(map(x->x.fit, bees))

    for i=1:No
        j = roulettSelect(bees, sum_f)
        updateBee!(bees[j], bees[getK(j, N)], f, bounds)
    end
end

function scoutPhase!(bees, f, genPos::Function, limit::Int)
    bees_scout = filter(x->x.t >= limit, bees)
    for i in 1:length(bees_scout)
        bees_scout[i].sol.x = genPos()
        bees_scout[i].sol.f = f(bees_scout[i].sol.x)
        bees_scout[i].t = 0
    end

    return length(bees_scout)
end

function getBestBee(bees)
    best = bees[1]
    for bee in bees
        if bee.sol.f < best.sol.f
            best = bee
        end
    end
 
    return best
end

function chooseBest(bees, best)
    bee_cand = getBestBee(bees) 
    if bee_cand.sol.f < best.sol.f
        return deepcopy(bee_cand)
    end

    return best
end

function initialbees(f, N, bounds)
     P = generate_population(f, N, bounds)

     return [ Bee(sol) for sol in P ]
end 

