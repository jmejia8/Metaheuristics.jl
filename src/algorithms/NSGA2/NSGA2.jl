include("crowding-distance.jl")

abstract type AbstractNSGA <: AbstractParameters end

# Abstracts for algorithm parameter
mutable struct NSGA2 <: AbstractNSGA
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
end

"""
    NSGA2(;
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

To use NSGA2, the output from the objective function should be a 3-touple
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

# define the parameters (use `NSGA2()` for using default parameters)
nsga2 = NSGA2(N = 100, p_cr = 0.85)

# optimize
status = optimize(f, bounds, nsga2)

# show results
display(status)
```

"""
function NSGA2(;
    N = 100,
    η_cr = 20,
    p_cr = 0.9,
    η_m = 20,
    p_m = -1,
    information = Information(),
    options = Options(),
)

    parameters = NSGA2(N, promote( Float64(η_cr), p_cr, η_m, p_m )...)
    Algorithm(
        parameters,
        information = information,
        options = options,
    )

end



function update_state!(
    status::State,
    parameters::AbstractNSGA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )


    I = randperm(parameters.N)
    J = randperm(parameters.N)
    Q = zeros(parameters.N, size(problem.bounds, 2))
    for i = 1:2:parameters.N

        pa = tournament_selection(status.population, I[i])
        pb = tournament_selection(status.population, J[i])

        c1, c2 = GA_reproduction(get_position(pa),
                                 get_position(pb),
                                 problem.bounds;
                                 η_cr = parameters.η_cr,
                                 p_cr = parameters.p_cr,
                                 η_m = parameters.η_m,
                                 p_m = parameters.p_m)
        Q[i,:] = c1
        Q[i+1,:] = c2
    end

    append!(status.population, create_solutions(Q, problem))

    # non-dominated sort, crowding distance, elitist removing
    environmental_selection!(status.population, parameters)

end

function reproduction(pa, pb, parameters::AbstractNSGA, problem)
    # crossover and mutation
    c1, c2 = GA_reproduction(get_position(pa),
                             get_position(pb),
                             problem.bounds;
                             η_cr = parameters.η_cr,
                             p_cr = parameters.p_cr,
                             η_m = parameters.η_m,
                             p_m = parameters.p_m)

    # evaluate offspring
    create_solution(c1, problem), create_solution(c2, problem) 
end

function environmental_selection!(population, parameters::AbstractNSGA)
    truncate_population!(population, parameters.N)
end

function initialize!(
    status,
    parameters::AbstractNSGA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    D = size(problem.bounds, 2)

    if parameters.p_m < 0.0
        parameters.p_m = 1.0 / D
    end

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    status = gen_initial_state(problem,parameters,information,options,status)

    status

end

function final_stage!(
    status::State,
    parameters::AbstractNSGA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()

end


function tournament_selection(P, a = rand(1:length(P)))
    # chose two different solutions at random
    b = rand(1:length(P))
    while a == b 
        b = rand(1:length(P))
    end

    # perform selection
    P[a].rank < P[b].rank || (P[a].rank == P[b].rank && P[a].crowding > P[b].crowding ) ? P[a] : P[b]
end

