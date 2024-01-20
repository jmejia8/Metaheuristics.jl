abstract type AbstractDifferentialEvolution <: AbstractParameters end

mutable struct DE <: AbstractDifferentialEvolution
    N::Int
    F::Float64
    CR::Float64
    CR_min::Float64
    CR_max::Float64
    F_min::Float64
    F_max::Float64
    strategy::Symbol
end

include("epsilonDE.jl")



"""
    DE(;
        N  = 0,
        F  = 1.0,
        CR = 0.5,
        strategy = :rand1,
        information = Information(),
        options = Options()
    )

Parameters for Differential Evolution (DE) algorithm: step-size `F`,`CR` controlls the binomial
crossover, `N` is the population size. The parameter `strategy` is related to the variation
operator (`:rand1`, `:rand2`, `:best1`, `:best2`, `:randToBest1`).

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE())
+=========== RESULT ==========+
  iteration: 1000
    minimum: 0
  minimizer: [0.0, 0.0, 0.0]
    f calls: 30000
 total time: 0.0437 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE(N=50, F=1.5, CR=0.8))
+=========== RESULT ==========+
  iteration: 600
    minimum: 8.68798e-25
  minimizer: [3.2777877981303293e-13, 3.7650459509488005e-13, -7.871487597385812e-13]
    f calls: 30000
 total time: 0.0319 s
+============================+
```

"""
function DE(;
        N::Int = 0,
        F = 0.7,
        CR = 0.5,
        CR_min = CR,
        CR_max = CR,
        F_min = F,
        F_max = F,
        strategy::Symbol = :rand1,
        kargs...
    )


    parameters =
    DE(N, promote(F, CR, CR_min, CR_max, F_min, F_max)..., strategy)

    Algorithm(parameters; kargs...)

end


function update_state!(
        status,
        parameters::AbstractDifferentialEvolution,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    population = status.population

    F = parameters.F
    CR = parameters.CR


    rng = options.rng
    # stepsize
    if parameters.F_min < parameters.F_max
        F = parameters.F_min + (parameters.F_max - parameters.F_min) * rand(rng)
    end

    if parameters.CR_min < parameters.CR_max
        CR =
        parameters.CR_min + (parameters.CR_max - parameters.CR_min) * rand(rng)
    end

    new_vectors = reproduction(status, parameters, problem)

    # evaluate solutions
    new_solutions = create_solutions(new_vectors, problem,ε=options.h_tol)
    append!(status.population, new_solutions)

    environmental_selection!(status.population, parameters)
    status.best_sol = get_best(status.population)
end


function environmental_selection(population, parameters::AbstractDifferentialEvolution)
    @assert length(population) == 2*parameters.N

    new_solutions = population[parameters.N+1:end]
    population = population[1:parameters.N]

    survivals = Int[]

    # select survival
    for (i, h) in enumerate(new_solutions)
        if is_better(h, population[i], parameters)
            push!(survivals, parameters.N + i)
        else
            push!(survivals, i)
        end
    end
    return survivals
end


function environmental_selection!(population, parameters::AbstractDifferentialEvolution)
    mask = environmental_selection(population, parameters)
    ignored = ones(Bool, length(population))
    ignored[mask] .= false
    deleteat!(population, ignored)

    return
end

function initialize!(
        status,
        parameters::AbstractDifferentialEvolution,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    D = getdim(problem)



    if parameters.N <= 5
        parameters.N = 10 * D
    end

    if parameters.CR < 0 || parameters.CR > 1
        parameters.CR = 0.5
        options.debug &&
        @warn("CR should be from interval [0,1]; set to default value 0.5")
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
        @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    return gen_initial_state(problem,parameters,information,options,status)
end

function final_stage!(
        status,
        parameters::AbstractDifferentialEvolution,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()
end


function reproduction(status, parameters::AbstractDifferentialEvolution, problem)
    @assert !isempty(status.population)

    N = parameters.N
    D = length(get_position(status.population[1]))

    strategy = parameters.strategy
    xBest = get_position(status.best_sol)
    population = status.population
    F = parameters.F
    CR = parameters.CR

    X = zeros(eltype(xBest), N,D)

    for i in 1:N
        x = get_position(population[i])
        u = DE_mutation(population, F, strategy, 1)
        v = DE_crossover(x, u, CR)
        evo_boundary_repairer!(v, xBest, problem.search_space)
        X[i,:] = _fix_type(v, problem.search_space)
    end

    X 
end

is_better(a, b, parameters::AbstractDifferentialEvolution) = is_better(a, b)
