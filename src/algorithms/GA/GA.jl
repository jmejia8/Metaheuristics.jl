abstract type AbstractGA <: AbstractParameters end

mutable struct GA{T1, T2, T3, T4, T5} <: AbstractGA
    initializer::T1
    selection::T2
    crossover::T3
    mutation::T4
    environmental_selection::T5
end

"""
    GA(;
        N = 100,
        p_mutation  = 1e-5,
        p_crossover = 0.5,
        initializer = RandomInBounds(),
        selection   = TournamentSelection(),
        crossover   = UniformCrossover(),
        mutation    = BitFlipMutation(),
        environmental_selection = ElitistReplacement()
        )

Just another Genetic Algorithm Framework.

### Parameters:

- `N` is the population size.
- `p_mutation`  mutation probability (gen).
- `p_crossover` crossover probability (gen).
- `initializer` parameters for the population initializer.
- `selection` parameters for the selection operator.
- `crossover` parameters for the crossover operators. 
- `mutation` parameters for the mutation operator.
- `environmental_selection` parameters for the replacement method.

### Example: Binary Encoding

```jldoctest
julia> f(x) = sum(x) / length(x)
f (generic function with 1 method)

julia> dim = 10;

julia> optimize(f, repeat([false, true], 1, dim), GA())
+=========== RESULT ==========+
  iteration: 500
    minimum: 0
  minimizer: Bool[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    f calls: 50000
 total time: 0.0523 s
stop reason: Maximum number of iterations exceeded.
+============================+
```

### Example: Permutation Encoding

```jldoctest
julia> f(x) = sum(abs.(x .- (length(x):-1:1.0)))
f (generic function with 1 method)

julia> perm_size = 10;

julia> optimize(f, zeros(Int,2,perm_size), ga)
+=========== RESULT ==========+
  iteration: 500
    minimum: 0
  minimizer: [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    f calls: 49900
 total time: 0.6230 s
stop reason: Maximum number of iterations exceeded.
+============================+
```

### Example: Integer Encoding

```jldoctest

```
"""
function GA(;
        N = 100,
        p_mutation  = 1e-5,
        p_crossover = 0.5,
        initializer = RandomInBounds(;N),
        selection   = TournamentSelection(;N),
        crossover   = UniformCrossover(p = p_crossover),
        mutation    = BitFlipMutation(p = p_mutation),
        environmental_selection = ElitistReplacement(),
        options     = Options(),
        information = Information()
    )

    parameters = GA(initializer,
                    selection,
                    crossover,
                    mutation,
                    environmental_selection
                   )

    Algorithm(
        parameters,
        information = information,
        options = options
    )
end


function initialize!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    @assert parameters.initializer.N > 0

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 2options.iterations * parameters.initializer.N
    end
    

    initialize!(status,parameters.initializer,problem,information,options,args...;kargs...)
end


function update_state!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    population = status.population

    # parent selection
    parent_mask = selection(population, parameters.selection)
    # variation operators
    Q = crossover(population[parent_mask], parameters.crossover)
    mutation!(Q, parameters.mutation)
    # evaluate offspring
    offsprings = create_solutions(Q, problem)
    cb = get_best(status.population)
    # select solutions (survivals are saved in population)
    environmental_selection!(population, offsprings, parameters.environmental_selection)

    if is_better(cb, status.best_sol)
        # save elite solution
        status.best_sol = cb
    end
    
end



function final_stage!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    status.final_time = time()
    return
end

function stop_criteria!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    return
end
