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
fvals(bees::Vector{<:Bee}) = fval.(bees)
is_feasible(bee::Bee) = is_feasible(bee.sol)


minimum(st::State{Bee{xf_indiv}}) = st.best_sol.sol.f
minimum(st::State{Bee{xfgh_indiv}}) = st.best_sol.sol.f
minimizer(st::State{Bee{xf_indiv}}) = st.best_sol.sol.x
minimizer(st::State{Bee{xfgh_indiv}}) = st.best_sol.sol.x

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

function updateBee!(bee, bee2, problem)
    D = length(bee.sol.x)
    ϕ = -1.0 + 2.0rand()

    v = ϕ*(bee.sol.x - bee2.sol.x)

    x_new = bee.sol.x + v
    replace_with_random_in_bounds!(x_new, problem.search_space)

    new_sol = create_solution(x_new, problem)

    if is_better(new_sol, bee.sol) #fx < bee.sol.f
        bee.sol = new_sol
        bee.fit = fit(new_sol.f)
        bee.t = 0
    else
        bee.t += 1
    end
end

function getK(i, N)
    j = rand(1:N)
    while j == i
        j = rand(1:N)
    end
    j
end

function employedPhase!(bees, problem, Ne)
    N = length(bees)
    for i in randperm(N)[1:Ne]
        updateBee!(bees[i], bees[getK(i, N)], problem)
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

function outlookerPhase!(bees, problem, No::Int)
    N = length(bees)
    sum_f = sum(map(x->x.fit, bees))

    for i=1:No
        j = roulettSelect(bees, sum_f)
        updateBee!(bees[j], bees[getK(j, N)], problem)
    end
end

function scoutPhase!(bees, problem, genPos::Function, limit::Int)
    bees_scout = filter(x->x.t >= limit, bees)
    for i in 1:length(bees_scout)
        bees_scout[i].sol = create_solution(genPos(), problem)
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

#=
function initialbees(N, problem)
     P = generate_population(N, problem)

     return [ Bee(sol) for sol in P ]
end 
=#
