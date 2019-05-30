mutable struct Bee
    sol
    fit::Float64
    t::Int
end

Bees = Array{Bee, 1}

Bee(solution) = Bee(solution, fit(solution.f), 0)


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

function updateBee!(bee, bee2, f)
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

function employedPhase!(bees, f, Ne)
    N = length(bees)
    S = view(bees, rand(1:N, Ne))
    for bee in S
        updateBee!(bee, bees[rand(1:N)], f)
    end
end

function roulettSelect(bees, sum_f)
    r = rand()
    rs = 0.0

    for bee in bees
        rs += bee.fit / sum_f
        if r <= rs
            return bee
        end
    end

    return bees[end]
end

function outlookerPhase!(bees, f, No::Int)
    N = length(bees)
    sum_f = sum(map(x->x.fit, bees))

    for i=1:No
        bee = roulettSelect(bees, sum_f)
        updateBee!(bee, bees[rand(1:N)], f)
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

function getBest(bees::Bees)
    best = bees[1]
    for bee in bees
        if bee.sol.f < best.sol.f
            best = bee
        end
    end
 
    return best
end

function chooseBest(bees, best::Bee)
    bee_cand = getBest(bees) 
    if bee_cand.sol.f < best.sol.f
        return deepcopy(bee_cand)
    end

    return best
end

function initialbees(f, N, bounds)
     P = initializePop(f, N, length(bounds[1,:]), bounds[1,:], bounds[2,:])

     return [ Bee(sol) for sol in P ]
end 

function ABC(
        fobj::Function,
        bounds;
        N = 50,
        limit=10,
        iters = Inf,
        max_evals = 10000*size(bounds, 2),
        Ne = div(N+1, 2),
        No = div(N+1, 2),
    )

    D = size(bounds, 2)
    @inline genPos(D=D, a=bounds[1,:], b = bounds[2,:]) = initializeSol(D, a, b)

    bees = initialbees(fobj, N, bounds)
    nevals = length(bees)

    best = deepcopy(getBest(bees))
    t = 0

    while nevals < max_evals && t < iters#i=1:iters
        t += 1

        employedPhase!(bees, fobj, Ne)        
        outlookerPhase!(bees, fobj, No)

        best = chooseBest(bees, best)

        nevals += Ne + No + scoutPhase!(bees, fobj, genPos, limit)

    end

    println(bees)

    return best.sol.x, best.sol.f
end

