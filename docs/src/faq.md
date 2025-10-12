# FAQ

Answers to Frequently Asked Questions:

## How to solve combinatorial problems?

This package was initially developed for numerical optimization, but recent updates
can handle combinatorial problems. The Genetic Algorithm framework ([`GA`](@ref)) can be used.

See the [N-Queens](@ref) tutorial for a detailed example.

## How to use different search spaces?

Metaheuristics.jl supports multiple search space types:

- **BoxConstrainedSpace** (or `boxconstraints`): For continuous optimization with bounds.
  ```julia
  bounds = boxconstraints(lb = zeros(10), ub = ones(10))
  ```

- **PermutationSpace**: For permutation-based problems (e.g., TSP, scheduling).
  ```julia
  space = PermutationSpace(10)  # permutations of size 10
  ```

- **BitArraySpace**: For binary optimization problems.
  ```julia
  space = BitArraySpace(20)  # 20 binary variables
  ```

- **MixedSpace**: For problems with mixed variable types.

See [Search Spaces](@ref) in the API reference for more details.

## How to handle constraints?

Constraints should be returned by your objective function. Use arrays for inequality (`g(x) ≤ 0`) and equality (`h(x) = 0`) constraints:

```julia
function f(x)
    fx = sum(x.^2)           # objective value
    gx = [sum(x) - 1]        # inequality: sum(x) ≤ 1
    hx = [x[1] - x[2]]       # equality: x[1] = x[2]
    return fx, gx, hx
end
```

For unconstrained problems, return `[0.0]` for both `gx` and `hx`.

## How to choose between algorithms?

Algorithm selection depends on your problem type:

- **Single-objective, unconstrained**: ECA, DE, PSO are good starting points
- **Single-objective, constrained**: εDE, ECA, SHADE
- **Multi-objective**: NSGA2, NSGA3, SPEA2, SMS_EMOA
- **Many-objective (>3 objectives)**: NSGA3
- **Combinatorial**: GA with appropriate operators, BRKGA
- **Large-scale**: CSO (supports large-scale problems)

See the [Algorithms Index](#) for a complete comparison table.

## How to set stopping criteria?

Use the `Options` struct to control when optimization stops:

```julia
options = Options(
    iterations = 1000,        # max iterations
    f_calls_limit = 50000,    # max function evaluations
    f_tol = 1e-6,            # stop when |f(x) - f_optimum| < f_tol
    time_limit = 60.0        # max time in seconds
)
```

Custom stopping criteria can be implemented by overriding `stop_criteria!`.

## How to perform parallel evaluation?

For batch evaluation of multiple solutions:

```julia
function f_parallel(X)  # X is N×D matrix
    N = size(X, 1)
    fx = zeros(N)
    Threads.@threads for i in 1:N
        fx[i] = my_expensive_function(X[i,:])
    end
    return fx
end

options = Options(parallel_evaluation=true)
optimize(f_parallel, bounds, ECA(options=options))
```

Start Julia with multiple threads: `julia -t 4`

See the [Parallelization](@ref) tutorial for more details.

## How to use custom initial solutions?

Use `set_user_solutions!` to provide starting solutions:

```julia
algorithm = ECA()
x0 = [0.5, 0.5, 0.5]  # single solution
set_user_solutions!(algorithm, x0, f)

# Or multiple solutions
X0 = rand(10, 3)  # 10 solutions
set_user_solutions!(algorithm, X0, f)

result = optimize(f, bounds, algorithm)
```

## How to visualize results?

Use Plots.jl to visualize results:

```julia
using Plots

# Get population positions
X = positions(result)
scatter(X[:,1], X[:,2])

# Plot convergence
f_calls, best_f = convergence(result)
plot(f_calls, best_f)

# For multi-objective: Pareto front
A = pareto_front(result)
scatter(A[:,1], A[:,2])
```

See the [Visualization](@ref) section for more examples.

## What's the difference between `optimize` and `optimize!`?

- **`optimize`**: Creates a new algorithm state and runs optimization. Use this most of the time.
- **`optimize!`**: Continues optimization from the current state. Useful for warm-starting or adaptive optimization.

```julia
result = optimize(f, bounds, ECA())    # fresh start
optimize!(result, f, bounds)            # continue from result
```

## How to perform multi-objective optimization?

Return a vector of objective values:

```julia
function f(x)
    fx = [sum(x.^2), sum((x .- 1).^2)]  # two objectives
    gx = [0.0]  # no constraints
    hx = [0.0]
    return fx, gx, hx
end

result = optimize(f, bounds, NSGA2())
pf = pareto_front(result)  # get Pareto front
```

## How to select the best solution from a Pareto front?

Use Multi-Criteria Decision Making (MCDM) methods:

```julia
# Define preferences (weights)
w = [0.5, 0.5]

# Use MCDM method
best_sol = best_alternative(result, w, CompromiseProgramming())

# Or use JMcDM methods
best_sol = best_alternative(result, w, TopsisMethod())
```

See [Multi-Criteria Decision-Making](@ref) for more details.


