# Abstracts for algorithm parameters
mutable struct SPEA2 <: AbstractNSGA
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    fitness::Vector{Float64}
end

"""
    SPEA2(;
        N = 100,
        η_cr = 20,
        p_cr = 0.9,
        η_m = 20,
        p_m = 1.0 / D,
        information = Information(),
        options = Options(),
    )

Parameters for the metaheuristic NSGA-II.

Parameters:

- `N` Population size.
- `η_cr`  η for the crossover.
- `p_cr` Crossover probability.
- `η_m`  η for the mutation operator.
- `p_m` Mutation probability (1/D for D-dimensional problem by default).

To use SPEA2, the output from the objective function should be a 3-touple
`(f::Vector, g::Vector, h::Vector)`, where `f` contains the objective functions,
`g` and `h` are inequality, equality constraints respectively.

A feasible solution is such that `g_i(x) ≤ 0 and h_j(x) = 0`.


```julia
using Metaheuristics

# Dimension
D = 2

# Objective function
f(x) = ( x, [sum(x.^2) - 1], [0.0] ) 

# bounds
bounds = [-1 -1;
           1  1.0
        ]

# define the parameters (use `SPEA2()` for using default parameters)
nsga2 = SPEA2(N = 100, p_cr = 0.85)

# optimize
status = optimize(f, bounds, nsga2)

# show results
display(status)
```

"""
function SPEA2(;
    N = 100,
    η_cr = 20,
    p_cr = 0.9,
    η_m = 20,
    p_m = -1,
    kargs...
)

    parameters = SPEA2(N, promote( Float64(η_cr), p_cr, η_m, p_m )...,[])
    Algorithm(parameters; kargs...)
end



function update_state!(
    status::State,
    parameters::SPEA2,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )

    if isempty(parameters.fitness)
        parameters.fitness = compute_fitness(status.population)
    end

    Q = reproduction(status, parameters, problem, options)

    append!(status.population, create_solutions(Q, problem))

    # non-dominated sort, crowding distance, elitist removing
    environmental_selection!(status.population, parameters)
end



function environmental_selection!(population, parameters::SPEA2)
    Distance = pairwise_distances(population)
    fitness = compute_fitness(population, Distance)
    N = parameters.N

    next = fitness .< 1
    K = count(next)
    if K < N
        rank = sortperm(fitness)
        next[rank[1:N]] .= true
    elseif K > N
        del  = truncation(population[next], K-N, Distance[next,next])
        temp = findall(next)
        next[temp[del]] .= false
    end

    deleteat!(population, .!next)
    deleteat!(fitness, .!next)
    parameters.fitness = fitness

    return
end

function truncation(population, K, distance = pairwise_distances(population))
    del = zeros(Bool, length(population))
    while count(del) < K
        remain   = findall(.!del)
        #temp     = sortperm(distance[remain,remain], dims = 2)
        # rank     = sortslicesperm(temp, dims=1)
        nn = argmin(distance[remain,remain])
        del[remain[nn.I[1]]] = true
    end

    return del
end

function compute_fitness(population, Distance = pairwise_distances(population))
    N = length(population)
    dominate = zeros(Bool,N, N)
    for i in 1:N
        for j in i+1:N
            k = compare(population[i], population[j])
            dominate[i,j] = k == 1
            dominate[j,i] = k == 2
        end
    end

    S = sum(dominate,dims=2)[:,1]
    R = [ sum(S[dominate[:,i]]) for i in 1:N]

    distance = sort(Distance,dims=2)


    D = 1 ./ (distance[:,floor(Int,sqrt(N))] .+ 2)
    return R + D
end


