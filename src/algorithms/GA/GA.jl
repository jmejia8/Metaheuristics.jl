abstract type AbstractGA <: AbstractParameters end

mutable struct GA{T1, T2, T3, T4, T5} <: AbstractGA
    initializer::T1
    selection::T2
    crossover::T3
    mutation::T4
    environmental_selection::T5
end

iscompatible(search_space::BitArraySpace, algorithm::GA) = true
iscompatible(search_space::PermutationSpace, algorithm::GA) = true

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

### Example: Binary Encoding (default)

```jldoctest
julia> f(x) = sum(x) / length(x)
f (generic function with 1 method)

julia> dim = 10;

julia> optimize(f, BitArraySpace(dim), GA())
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

julia> ga = GA(;initializer = RandomPermutation(N=100), crossover=OrderCrossover(), mutation=SlightMutation());

julia> optimize(f, PermutationSpace(perm_size), ga)
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
julia> c = randn(10);

julia> A = randn(10, 10);

julia> f(x) = c'*x, A*x, [0.0]
f (generic function with 1 method)

julia> bounds = repeat([0,100], 1, 10)
2×10 Matrix{Int64}:
   0    0    0    0    0    0    0    0    0    0
 100  100  100  100  100  100  100  100  100  100

julia> ga = GA(;crossover=SBX(;bounds),
                mutation=PolynomialMutation(;bounds),
                environmental_selection=GenerationalReplacement()
            );

julia> result = optimize(f, bounds, ga)
+=========== RESULT ==========+
  iteration: 500
    minimum: -76.9031
  minimizer: [0, 10, 49, 100, 24, 0, 0, 0, 67, 76]
    f calls: 50000
  feasibles: 95 / 100 in final population
 total time: 0.6028 s
stop reason: Maximum number of iterations exceeded.
+============================+
```

### Example: Real Encoding

```jldoctest
julia> f, bounds, solutions = Metaheuristics.TestProblems.rastrigin();

julia> ga = GA(;mutation=PolynomialMutation(;bounds),
                crossover=SBX(;bounds),
                environmental_selection=GenerationalReplacement()
               );

julia> optimize(f, bounds, ga)
+=========== RESULT ==========+
  iteration: 500
    minimum: 5.91136e-09
  minimizer: [-5.446073166732363e-6, -1.7900876625850504e-7, 2.8548505431323723e-8, 2.1130514595980084e-8, -1.632381298278562e-8, 2.8216219650016283e-9, -3.2114953600427333e-7, 2.4222648522114125e-9, -5.5236545829928716e-9, -2.3479408274628334e-9]
    f calls: 49900
 total time: 0.5775 s
stop reason: Maximum number of iterations exceeded.
+============================+
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

# this method is called in src/optimize.jl
function get_parameters(
        f,
        search_space::BoxConstrainedSpace, # real/integer encoding
        ::Type{T}
    ) where T <: GA 
    
    N = 100
    GA(;
        N = N,
        initializer = RandomInBounds(;N),
        selection   = TournamentSelection(;N),
        crossover   = SBX(;bounds=search_space),
        mutation    = PolynomialMutation(;bounds=search_space),
        environmental_selection = ElitistReplacement()
      )
end


function get_parameters(
        f,
        search_space::PermutationSpace, # permutation-based encoding
        ::Type{T}
    ) where T <: GA

    N = 100
    GA(;
       N = N,
       initializer = RandomPermutation(;N),
       selection   = TournamentSelection(;N),
       crossover   = OrderCrossover(),
       mutation    = SlightMutation(),
       environmental_selection = ElitistReplacement()
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

