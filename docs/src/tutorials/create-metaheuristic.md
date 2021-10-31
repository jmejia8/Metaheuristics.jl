# Create Your Own Metaheuristic

## Introduction 

Firstly, you need to know what occurs when the [`optimize`](@ref) function is called.

### Optimization Process


1. **Initialization**: `status = initialize!(status, parameters, problem, information, options)`
   this function should initialize a `State` with population members according
   to the parameters provided.
2. **Main optimization loop**: while `status.stop == false` do
    - update population, parameters via `update_state!(status, parameters, problem, information, options)`,
    - and `stop_criteria!(status, parameters, problem, information, options)` will change `status.stop`.
3. **Final Stage**: When the loop in step 2 breaks, then a final function is called `final_stage!`
   for the final update of the state, e.g., delete infeasible solutions in population,
   get non-dominated solutions, etc. 


**Initialization**:

```julia
function initialize!(
                status, # an initiliazed State (if apply)
                parameters::AbstractParameters,
                problem,
                information,
                options,
                args...;
                kargs...
        )

    # initialize the stuff here
    return State(0.0, zeros(0)) # replace this
end
```

**Optimization Process**: In this step, the [`State`](@ref) is updated using the following
function which is called at each iteration/generation.

```julia
function update_state!(
        status,
        parameters::AbstractParameters,
        problem,
        information,
        options,
        args...;
        kargs...
)
    # update any element in State 
    return
end
```


**Final Step:**

```julia
function final_stage!(
        status,
        parameters::AbstractParameters,
        problem,
        information,
        options,
        args...;
        kargs...
)
    return
end
```

### The Algorithm Parameters

Any proposed algorithm, let's say "XYZ", uses different parameters, then it is suggested to store them in a
structure, e.g.:

```julia
# structure with algorithm parameters
mutable struct XYZ <: AbstractParameters
    N::Int # population size
    p_crossover::Float64 # crossover probability
    p_mutation::Float64 # mutation probability
end

# a "constructor" 
function XYZ(;N = 0, p_crossover = 0.9, p_mutation = 0.1)
    parameters = XYZ(N, p_crossover, p_mutation)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end
```


If you want to implement an algorithm outside of the `Metaheuristics` module, then
include explicitly the methods you require (or use the `Metaheuristics.` prefix)
as in Step 0, otherwise go to Step 1.


## Implementing a Simple Genetic Algorithm

The following steps describe how to implement a simple Genetic Algorithm.

### Step 0

Including stuff from `Metaheuristics` we need.

```julia
# base methods
using Metaheuristics
import Metaheuristics: initialize!, update_state!, final_stage!
import Metaheuristics: AbstractParameters, gen_initial_state, Algorithm, get_position
# genetic operators
import Metaheuristics: SBX_crossover, polynomial_mutation!, create_solution, is_better
import Metaheuristics: reset_to_violated_bounds!
```

### Step 1: The Parameters

Due to we are creating a simple Genetic Algorithm (GA), let's define the parameters for the GA.

```julia
# structure with algorithm parameters
mutable struct MyGeneticAlgorithm <: AbstractParameters
    N::Int # population size
    p_crossover::Float64 # crossover probability
    p_mutation::Float64 # mutation probability
end
```

```julia
function MyGeneticAlgorithm(;N = 100,
                            p_crossover = 0.9,
                            p_mutation = 0.1,
                            information = Information(),
                            options = Options()
    )
    parameters = MyGeneticAlgorithm(N, p_crossover, p_mutation)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end
```


### Step 2: Initialization

Initialize population, parameters and settings before the optimization process begins.
The most common initialization method is generating uniformly distribution random 
number in provided bounds. Here, [`gen_initial_state`](@ref) for that purpose. Note
that [`gen_initial_state`](@ref) require that `parameters.N` is defined.

```julia
function initialize!(
        status,
        parameters::MyGeneticAlgorithm,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    # gen_initial_state require that `parameters.N` is defined.
    return gen_initial_state(problem,parameters,information,options,status)

end
```

### Step 3: Evolve Population

Now, it is time to update (evolve) your population by using genetic operators: selection,
crossover, mutation and environmental selection.

```julia
function update_state!(
        status,
        parameters::MyGeneticAlgorithm,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    population = status.population
    N = parameters.N

    for i in 1:N
        # selection
        parent_1 = get_position(rand(population))
        parent_2 = get_position(rand(population))

        # generate offspring  via SBX crossover
        c,_ = SBX_crossover(parent_1, parent_2, problem.bounds, 20, parameters.p_crossover)

        # Mutate solution
        polynomial_mutation!(c, problem.bounds, 15, parameters.p_mutation)
        # Fix solution if necessary
        reset_to_violated_bounds!(c, problem.bounds)

        # crate the solution and evaluate fitness (x, f(x))
        offspring = create_solution(c, problem)

        push!(population, offspring)
    end

    # environmental selection
    sort!(population, lt = is_better, alg=PartialQuickSort(N))
    deleteat!(population, N+1:length(population))

end
```

### Step 4: After Evolution

This step is optional, but here is used to get the elite solution aka the best solution
found by our GA.

```julia
function final_stage!(
        status,
        parameters::MyGeneticAlgorithm,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    # first solution is the best one since population is ordered in previous step
    status.best_sol = status.population[1]
    status.final_time=time()
    return
end
```

### Step 5: Time to Optimize

Now, we are able to solve and optimization problem using our genetic algorithm.


!!! compat "Optimization Problems"
    As you can see, `MyGeneticAlgorithm` was not restricted to any kind of optization problems,
    however works for constrained, unconstrained single- and multi-objective problems; why?
    The method `gen_initial_state` creates a [`State`](@ref) according to the output
    of the objective function `f`, whilst `is_better` is comparing solutions according
    to the solution type.

```julia
function main()
    # test problem
    f, bounds, _ = Metaheuristics.TestProblems.rastrigin()

    # optimize and get the results
    res = optimize(f, bounds, MyGeneticAlgorithm())
    display(res)
end

main()
```

**Output:**

```
+=========== RESULT ==========+
  iteration: 500
    minimum: 1.41152e-06
  minimizer: [6.98505513305995e-6, -1.1651666613615994e-5, -4.343967193003195e-6, 3.567365134464557e-5, 1.3393840640183734e-5, 5.591802709915942e-5, -1.477407456986382e-5, 6.325103756718973e-6, 1.9153467328726614e-5, 4.132106648380982e-5]
    f calls: 50000
 total time: 1.3685 s
+============================+
```

See [`optimize`](@ref) for more information.

### Exercises

1. Test your algorithm on a multi-objective optimization problem. Suggestion: change `rastrigin`
   by `ZDT1`.
2. Implement an interest metaheuristic and make a PR to the [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl) on github.
