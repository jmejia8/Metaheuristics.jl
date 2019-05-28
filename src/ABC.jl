using Printf
include("structures.jl")
include("operators.jl")
include("tools.jl")

mutable struct Bee
    sol
    fitness :: AbstractFloat
    t :: Integer
end

Bee(data) = Bee(data, 0.0, 0)


function Base.copy(s::Bee)
    return Bee(copy(s.sol.x), s.fitness, s.t)
end


@inline function update_fitness!(bee::Bee, f)
    bee.fitness = fitness(f(bee.sol.x))
end

@inline function update_fitness!(Bees::Array, f)
    for bee in Bees
        update_fitness!(bee, f)
    end
end

function fitness(val)
    if val >= 0.0
        return 1.0/(1.0+val)
    end
    
    return 1.0+abs(val)
end

function update_bee!(bee, beem::Bee, f)
    D = length(bee.sol.x)
    ϕ = -1.0 + 2.0rand()
    j = rand(1:D)

    x_new = ϕ*(bee.sol.x[j] - beem.sol.x[j])

    bee.sol.x[j] += x_new
    fit_new = fitness(f(bee.sol.x))

    if fit_new > bee.fitness
        bee.fitness = fit_new
        bee.t = 0
    else
        bee.sol.x[j] -= x_new
        bee.t += 1
    end
end

function update_employed!(Bees, f)
    N = length(Bees)
    for bee in Bees
        update_fitness!(bee, f)
        update_bee!(bee, Bees[rand(1:N)], f)
    end
end

function roulette_select(Bees)
    sf = sum(map(x->x.fitness, Bees))
    r = rand()
    rs = 0.0

    for bee in Bees
        rs += bee.fitness / sf
        if r <= rs
            return bee
        end
    end

    return Bees[end]
end

function update_outlook!(Bees, f, No::Integer)
    N = length(Bees)
    for i=1:No
        bee = roulette_select(Bees)
        update_fitness!(bee, f)
        update_bee!(bee, Bees[rand(1:N)], f)
    end
end

function update_scout!(Bees, f, genBee::Function, limit::Integer)
    bees_scout = filter(x->x.t >= limit, Bees)
    for bee in bees_scout
        bee.sol.x = genBee()
        bee.t = 0
    end
    update_fitness!(bees_scout, f)
end

function update_scout!(Bees, f, genBee::Function, limit::Integer)
    bees_scout = filter(x->x.t >= limit, Bees)
    for bee in bees_scout
        bee.sol.x = genBee()
        bee.t = 0
    end
end

function getBest(Bees)
    best = Bees[1]
    for bee in Bees
        if bee.fitness > best.fitness
            best = bee
        end
    end
    return best
end

function updateBest!(Bees, bee_best::Bee)
    bee_cand = getBest(Bees)
    if bee_best.fitness < bee_cand.fitness
        bee_best = deepcopy(bee_cand)
    end
end

# function updateBest!(Bees, bee_best::Bee)
#     bee_best = deepcopy(bee_cand)
# end

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
        No =10
    )

    D = size(bounds, 2)
    @inline genBee(D=D, a=bounds[1,:], b = bounds[2,:]) = initializeSol(D, a, b)

    Bees = initialBees(fobj, N, bounds)

    best = getBest(Bees)
    for i=1:iters
        update_employed!(Bees, fobj)
        update_outlook!(Bees, fobj, No)
        updateBest!(Bees, best)
        update_scout!(Bees, fobj, genBee, limit)
    end
    return best.sol.x, best.sol.f
end

f(x) = sum(x.^2)

ABC(f, [-10.0 -10 -10; 10 10 10])