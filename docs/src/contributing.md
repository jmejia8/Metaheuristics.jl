# Contributing

Please, be free to send me your PR, issue or any comment about this package for Julia.

### Algorithm Structure

When you call the [`optimize`](@ref) function, the following steps are carried out:

1. Initialization: `initialize!(problem, engine, parameters, status, information, options)`
2. Main optimization loop: while `status.stop == false` do
    - update population, parameters via `update_state!(problem, engine, parameters, status, information, options, iteration)`, and 
    - mainly set `status.stop = engine.stop_criteria(status, information, options)`
3. When the loop in step 2 beaks, then a final function is called `final_stage!` in order
   to update or refine the final state, e.g., delete infeasible solutions in population,
   get non-dominated solutions, etc. 

**Initialization**:

```julia
function initialize!(
    problem::Problem,
    engine::Engine,
    parameters::Any,
    status::State,
    information::Information,
    options::Options,
   )
    # initialize parameters, population, etc.
end
```

**Optimization Process**: In this step, the [`State`](@ref) is updated using the following
function which is called at each iteration/generation.

```julia
function update_state!(
        problem::Problem,
        engine::Engine,
        parameters::Any,
        status::State,
        information::Information,
        options::Options,
        iteration::Int,
       )
    # update any element in State 
end
```


**Final Step:**

```julia
function final_stage!(status::State, information::Information, options::Options)
    # used to a final update of the status. 
end
```

### Parameter

Any proposed algorithm, let's say "XYZ", uses different parameters, then it is suggested to store them in a
structure, e.g.:

```julia

# structure
mutable struct XYZ <: AbstractAlgorithm
    N::Int # population size
    p_crossover::Float64 # crossover probability
    p_mutation::Float64 # mutation probability
end

# constructor like
function XYZ(;N = 0, p_crossover = 0.9, p_mutation = 0.1)
    parameters = XYZ(N, p_crossover, p_mutation)

    Algorithm(
        parameters,
        initialize! = initialize!,
        update_state! = update_state!,
        is_better = is_better,
        stop_criteria = stop_check,
        final_stage! = final_stage!,
        information = information,
        options = options,
    )
end
```

