using Printf
include("structures.jl")
include("operators.jl")
include("tools.jl")

mutable struct Bee
    sol
    fit::Float64
    t::Int
end

Bee(solution) = Bee(solution, 0.0, 0)


function Base.copy(s::Bee)
    return Bee(copy(s.sol.x), s.fit, s.t)
end


@inline function update_fit!(bee::Bee)
    bee.fit = fit(bee.sol.f)
end

@inline function update_fit!(Bees::Array)
    for bee in Bees
        update_fit!(bee)
    end
end

function fit(fx)
    if fx >= 0.0
        return 1.0/(1.0+fx)
    end
    
    1.0 + abs(fx)
end

function updateBee!(bee, beem::Bee, f)
    D = length(bee.sol.x)
    ϕ = -1.0 + 2.0rand()

    v = ϕ*(bee.sol.x - beem.sol.x)

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

function employedPhase!(Bees, f)
    N = length(Bees)
    for bee in Bees
        updateBee!(bee, Bees[rand(1:N)], f)
    end
end

function roulettSelect(Bees, sum_f)
    r = rand()
    rs = 0.0

    for bee in Bees
        rs += bee.fit / sum_f
        if r <= rs
            return bee
        end
    end

    return Bees[end]
end

function outlookerPhase!(Bees, f, No::Integer)
    N = length(Bees)
    sum_f = sum(map(x->x.fit, Bees))

    for i=1:No
        bee = roulettSelect(Bees, sum_f)
        updateBee!(bee, Bees[rand(1:N)], f)
    end
end

function scoutPhase!(Bees, f, genBee::Function, limit::Integer)
    bees_scout = filter(x->x.t >= limit, Bees)
    for bee in bees_scout
        bee.sol.x = genBee()
        bee.t = 0
    end
end

function getBest(Bees::Array{Bee})
    best = Bees[1]
    for bee in Bees
        if bee.sol.f < best.sol.f
            best = bee
        end
    end
 
    return best
end

function chooseBest(Bees, best::Bee)
    bee_cand = getBest(Bees) 
    if bee_cand.sol.f < best.sol.f
        return deepcopy(bee_cand)
    end

    return best
end

function initialBees(f, N, bounds)
     P = initializePop(f, N, length(bounds[1,:]), bounds[1,:], bounds[2,:])

     return [ Bee(sol) for sol in P ]
end 

function ABC(
        fobj::Function,
        bounds;
        N = 100,
        limit=10,
        iters = 1000,
        Ne = div(N+1, 2),
        No = div(N+1, 2),
    )

    D = size(bounds, 2)
    @inline genBee(D=D, a=bounds[1,:], b = bounds[2,:]) = initializeSol(D, a, b)

    Bees = initialBees(fobj, N, bounds)

    best = deepcopy(getBest(Bees))
    
    for i=1:iters
        employedPhase!(Bees, fobj)
        outlookerPhase!(Bees, fobj, No)
        best = chooseBest(Bees, best)
        scoutPhase!(Bees, fobj, genBee, limit)
    end
    return best.sol.x, best.sol.f
end

f(x) = sum(x.^2)

@time ABC(f, [-10.0 -10 -10; 10 10 10])