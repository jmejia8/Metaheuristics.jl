mutable struct SMS_EMOA <: AbstractNSGA
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    n_samples::Int
end


include("update-population.jl")

"""
    SMS_EMOA(;
        N = 100,
        η_cr = 20,
        p_cr = 0.9,
        η_m = 20,
        p_m = 1.0 / D,
        n_samples = 10_000,
        information = Information(),
        options = Options(),
    )

Parameters for the metaheuristic SMS-EMOA.

Parameters:

- `N` Population size.
- `η_cr`  η for the crossover.
- `p_cr` Crossover probability.
- `η_m`  η for the mutation operator.
- `p_m` Mutation probability (1/D for D-dimensional problem by default).
- `n_samples` number of samples to approximate hypervolume in many-objective (M > 2).

To use SMS_EMOA, the output from the objective function should be a 3-touple
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

# define the parameters (use `SMS_EMOA()` for using default parameters)
sms_emoa = SMS_EMOA(N = 100, p_cr = 0.85)

# optimize
status = optimize(f, bounds, sms_emoa)

# show results
display(status)
```

"""
function SMS_EMOA(;
    N = 100,
    η_cr = 20,
    p_cr = 0.9,
    η_m = 20,
    p_m = -1,
    n_samples = 10_000,
    kargs...
)

    parameters = SMS_EMOA(N, promote( Float64(η_cr), p_cr, η_m, p_m)..., n_samples)
    Algorithm(parameters; kargs...)

end



function update_state!(
    status::State,
    parameters::SMS_EMOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )


    I = randperm(parameters.N)
    J = randperm(parameters.N)

    if options.parallel_evaluation
        Q = zeros(parameters.N, getdim(problem))
    end 

    for i = 1:parameters.N
        pa = tournament_selection(status.population, options, I[i])
        pb = tournament_selection(status.population, options, J[i])

        # crossover
        _, c = SBX_crossover(get_position(pa),
                             get_position(pb),
                             problem.search_space,
                             parameters.η_cr,
                             parameters.p_cr)
       
        # mutation
        polynomial_mutation!(c, problem.search_space, parameters.η_m, parameters.p_m)
       
        # repair solutions if necesary
        reset_to_violated_bounds!(c, problem.search_space)

        if options.parallel_evaluation
            Q[i,:] = c
            continue
        end

        # evaluate offspring
        offspring = create_solution(c, problem)
       
        # Reduce: save offspring in population and reduce population according to
        # front contribution
        update_population!(status.population, offspring, parameters.n_samples)
    end

    if !options.parallel_evaluation
        return
    end
    append!(status.population, create_solutions(Q, problem))
    environmental_selection!(status.population, parameters)
end


function environmental_selection!(population, parameters::SMS_EMOA)
    if length(population) <= parameters.N
        return
    end

    new_solutions = population[parameters.N+1:end]
    population = population[1:parameters.N]

    for (i, offspring) in enumerate(new_solutions)
        # Reduce: save offspring in population and reduce population according to
        # front contribution
        update_population!(population, offspring, parameters.n_samples)
    end
end

function initialize!(
    status,
    parameters::SMS_EMOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...)

    D = getdim(problem)

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
    parameters::SMS_EMOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()

    #status.best_sol = get_pareto_front(status.population, is_better)

end

