mutable struct Bee
    data :: Vector
    fitness :: AbstractFloat
    count :: Integer
end

function Bee(data::Vector)
    Bee(data, zero(eltype(data)), zero(Integer))
end


function Base.copy(s::Bee)
    return Bee(copy(s.data), s.fitness, s.count)
end


@inline function update_fitness!(bee::Bee, g::Function)
    bee.fitness = fitness(g(bee.data))
end

@inline function update_fitness!(bees::Bees, g::Function)
    for bee in bees
        update_fitness!(bee, g)
    end
end

function update_bee!(bee::Bee, beem::Bee, g::Function)
    dim = length(bee.data)
    phi = 2.0*(rand()-0.5)
    j = rand(1:dim)

    val = phi*(bee.data[j] - beem.data[j])

    bee.data[j] += val
    fitnew = fitness(g(bee.data))

    if fitnew > bee.fitness
        bee.fitness = fitnew
        bee.count = 0
    else
        bee.data[j] -= val
        bee.count += 1
    end
end

function update_employed!(bees::Bees, g::Function)
    N = length(bees)
    for bee in bees
        update_bee!(bee, bees[rand(1:N)], g)
    end
end

function update_employed!(bees::Bees, g::Function)
    N = length(bees)
    for bee in bees
        update_fitness!(bee, g)
        update_bee!(bee, bees[rand(1:N)], g)
    end
end

function roulette_select(bees::Bees)
    sf = sum(map(x->x.fitness, bees))
    r = rand()
    rs = zero(r)

    for bee in bees
        rs += bee.fitness/sf
        if r<=rs
            return bee
        end
    end
    return bees[end]
end

function update_outlook!(bees::Bees, g::Function, No::Integer)
    N = length(bees)
    for i=1:No
        bee = roulette_select(bees)
        update_bee!(bee, bees[rand(1:N)], g)
    end
end

function update_outlook!(bees::Bees, g::Function, No::Integer)
    N = length(bees)
    for i=1:No
        bee = roulette_select(bees)
        update_fitness!(bee, g)
        update_bee!(bee, bees[rand(1:N)], g)
    end
end

function update_scout!(bees::Bees, g::Function, init::Function, limit::Integer)
    bees_scout = filter(x->x.count >= limit, bees)
    for bee in bees_scout
        bee.data = init()
        bee.count = 0
    end
    update_fitness!(bees_scout, g)
end

function update_scout!(bees::Bees, g::Function, init::Function, limit::Integer)
    bees_scout = filter(x->x.count >= limit, bees)
    for bee in bees_scout
        bee.data = init()
        bee.count = 0
    end
end

function find_best(bees::Bees)
    best = bees[1]
    for bee in bees
        if bee.fitness > best.fitness
            best = bee
        end
    end
    return best
end

function update_best!(bees::Bees, bee_best::Bee)
    bee_cand = find_best(bees)
    if bee_best.fitness < bee_cand.fitness
        copyto!(bee_best, bee_cand)
    end
end

function update_best!(bees::Bees, bee_best::Bee)
    copyto!(bee_best, find_best(bees))
end

function ABC(
        fobj::Function,
        bounds;
        N = 100,
        limit=10,
        iters = 1000,
        No
    )

    D = size(bounds, 2)
    Population = initializePop(fobj, N, D, bounds[1,:], bounds[2,:])

    for i=1:iters
        update_employed!(Population, fobj)
        update_outlook!(Population, fobj, No)
        update_best!(Population, best)
        update_scout!(Population, fobj, abc.init, limit)
    end
    return abc.best.data
end

function fitness(val)
    if val >= 0.0
        return 1.0/(1.0+val)
    end
    
    return 1.0+abs(val)
end

