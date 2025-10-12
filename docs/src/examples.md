# Examples

This package provides different tools for optimization. Hence, this section gives different
examples for using the implemented `Metaheuristics`.

## Single-Objective Optimization


Firstly import this package

```@example SingleObjective
using Metaheuristics
```
Now, let us define the objective function to be minimized:

```@example SingleObjective
f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )
```

The search space (a.k.a. box-constraints) can be defined as follows:

```@example SingleObjective
bounds = boxconstraints(lb = -5ones(10), ub = 5ones(10))
```

!!! note "Search Space Definition"
    `boxconstraints` is an alias for `BoxConstrainedSpace`. The matrix format
    `bounds = [-5ones(10) 5ones(10)]'` is also supported but using named parameters
    is recommended for clarity.

It is possible to provide some information on the minimization problem.
Let's provide the true optimum to stop the optimizer when a tolerance `f_tol` is satisfied.

```@example SingleObjective
information = Information(f_optimum = 0.0)
```

Generic options or settings (e.g. budget limitation, tolerances, etc) can be provided as follows:

```@example SingleObjective
options = Options(f_calls_limit = 9000*10, f_tol = 1e-5, seed=1)
```

Now, we can provide the Information and Options to the optimizer (ECA in this example).

```@example SingleObjective
algorithm = ECA(information = information, options = options)
```

Now, the optimization is performed as follows:

```@example SingleObjective
result = optimize(f, bounds, algorithm)
```

The minimum and minimizer:

```@example SingleObjective
minimum(result)
```

```@example SingleObjective
minimizer(result)
```


!!! compat "Second run is faster in Julia"
    As you may know, the second run can be faster.


## Constrained Optimization

It is common that optimization models include constraints that must be satisfied. For example:
[The Rosenbrock function constrained to a disk](https://en.wikipedia.org/wiki/Test_functions_for_optimization)

Minimize:

```math
{\displaystyle f(x,y)=(1-x)^{2}+100(y-x^{2})^{2}}
```
subject to:

```math
{\displaystyle x^{2}+y^{2}\leq 2}
```
where $-2 \leq x,y \leq 2$.

In `Metaheuristics.jl`, a feasible solution is such that $g(x) \leq 0$ and $h(x) \approx 0$.
Hence, in this example the constraint is given by $g(x) = x^2 + y^2 - 2 \leq 0$.
Moreover, the equality and inequality constraints must be saved into  `Array`s.

!!! compat "Constraints handling"
    In this package, if the algorithm was not designed for constrained optimization,
    then solutions with the lower constraint violation sum will be preferred.

```@repl
using Metaheuristics

function f(x)
    x,y = x[1], x[2]

    fx = (1-x)^2+100(y-x^2)^2
    gx = [x^2 + y^2 - 2] # inequality constraints
    hx = [0.0] # equality constraints

    # order is important
    return fx, gx, hx
end

bounds = boxconstraints(lb = [-2.0, -2.0], ub = [2.0, 2.0])

optimize(f, bounds, ECA(N=30, K=3))
```

## Multiobjective Optimization

To implement a multiobjective optimization problem and solve it, you can proceed as usual. Here,
you need to provide constraints if they exist, otherwise put `gx = [0.0]; hx = [0.0];`
to indicate an unconstrained multiobjective problem.

```@repl
using UnicodePlots # to visualize in console (optional)
using Metaheuristics

function f(x)
    # objective functions
    v = 1.0 + sum(x .^ 2)
    fx1 = x[1] * v
    fx2 = (1 - sqrt(x[1])) * v

    fx = [fx1, fx2]

    # constraints
    gx = [0.0] # inequality constraints
    hx = [0.0] # equality constraints

    # order is important
    return fx, gx, hx
end

bounds = boxconstraints(lb = zeros(30), ub = ones(30))

optimize(f, bounds, NSGA2())
```


## Bilevel Optimization

Bilevel optimization problems can be solved by using the package
[BilevelHeuristics.jl](https://github.com/jmejia8/BilevelHeuristics.jl) which extends
 `Metaheuristics.jl` for handling those hierarchical problems.

Defining objective functions corresponding to the BO problem.

**Upper level (leader problem):**

```julia
using BilevelHeuristics

F(x, y) = sum(x.^2) + sum(y.^2)
bounds_ul = boxconstraints(lb = -ones(5), ub = ones(5)) 
```

**Lower level (follower problem):**

```julia
f(x, y) = sum((x - y).^2) + y[1]^2
bounds_ll = boxconstraints(lb = -ones(5), ub = ones(5))
```
**Approximate solution:**

```julia
res = optimize(F, f, bounds_ul, bounds_ll, BCA())
```

**Output:**
```
+=========== RESULT ==========+
  iteration: 108
    minimum: 
          F: 4.03387e-10
          f: 2.94824e-10
  minimizer: 
          x: [-1.1460768817533927e-5, 7.231706879604178e-6, 3.818596951258517e-6, 2.294324313691869e-6, 1.8770952450067828e-6]
          y: [1.998748659975197e-6, 9.479307908087866e-6, 6.180041276047425e-6, -7.642051857319683e-6, 2.434166021682429e-6]
    F calls: 2503
    f calls: 5062617
    Message: Stopped due UL function evaluations limitations. 
 total time: 26.8142 s
+============================+
```

See [BilevelHeuristics](https://jmejia8.github.io/BilevelHeuristics.jl/dev/) documentation
for more information.


## Decision-Making

Although Metaheuristics is focused on the optimization part, some decision-making algorithms
are available in this package (see [Multi-Criteria Decision-Making](@ref)).

The following example shows how to perform *a posteriori* decision-making.

```julia-repl
julia> # load the problem
julia> f, bounds, pf = Metaheuristics.TestProblems.ZDT1();

julia> # perform multi-objective optimization
julia> res = optimize(f, bounds, NSGA2());

julia> # user preferences
julia> w = [0.5, 0.5];

julia> # set the decision-making algorithm
julia> dm_method = CompromiseProgramming(Tchebysheff())

julia> # find the best decision
julia> sol = best_alternative(res, w, dm_method)
(f = [0.38493217206706115, 0.38037042164979956], g = [0.0], h = [0.0], x = [3.849e-01, 7.731e-06, …, 2.362e-07])
```

## Providing Initial Solutions

Sometimes you may need to use the starter solutions you need before the optimization
process begins, well, this example illustrates how to do it.

```@repl
using Metaheuristics # hide
f(x) = abs(x[1]) + x[2]  + x[3]^2 # objective function
algo  = ECA(N = 61); # optimizer

# one solution can be provided
x0 = [0.5, 0.5, 0.5];

set_user_solutions!(algo, x0, f);

# or multiple solutions can be given
X0 = rand(30, 3); # 30 solutions with dim 3

set_user_solutions!(algo, X0, f);
optimize(f, boxconstraints(lb = [0, 0, 0], ub = [1, 1, 1]), algo)
```


## Batch Evaluation

Evaluating multiple solutions at the same time can reduce computational time. To do that,
define your function on an input `N x D` matrix and function values into matrices with outcomes
in rows for all `N` solutions. Also, you need to put `parallel_evaluation=true` in the [`Options`](@ref)
to indicate that your `f` is prepared for parallel evaluations.

```julia
f(X) = begin
    fx = sum(X.^2, dims=2)       # objective function ∑x²
    gx = sum(X.^2, dims=2) .-0.5 # inequality constraints ∑x² ≤ 0.5
    hx = zeros(0,0)              # equality constraints
    fx, gx, hx
end

options = Options(parallel_evaluation=true)

res = optimize(f, boxconstraints(lb = -10ones(5), ub = 10ones(5)), ECA(options=options))
```

See [Parallelization](@ref) tutorial for more details.


## Modifying an Existing Metaheuristic

You may need to modify one of the implemented metaheuristics to improve the
algorithm performance or test new mechanisms. This example illustrates how to do it.


!!! warning "Modifying algorithms could break stuff"
    Be cautious when modifying a metaheuristic due to those changes will overwrite the default
    method for that metaheuristic.


Let's assume that we want to modify the stop criteria for `ECA`. See [Contributing](@ref) 
for more details.

```julia
using Metaheuristics
import LinearAlgebra: norm

# overwrite method
function Metaheuristics.stop_criteria!(
        status,
        parameters::ECA, # It is important to indicate the modified Metaheuristic 
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if status.stop
        # nothing to do
        return
    end

    # Diversity-based stop criteria

    x_mean = zeros(length(status.population[1].x))
    for sol in status.population
        x_mean += sol.x
    end
    x_mean /= length(status.population)
    
    distances_mean = sum(sol -> norm( x_mean - sol.x ), status.population)
    distances_mean /= length(status.population)

    # stop when solutions are close enough to the geometrical center
    new_stop_condition = distances_mean <= 1e-3

    status.stop = new_stop_condition

    # (optional and not recommended) print when this criterium is met
    if status.stop
        @info "Diversity-based stop criterium"
        @show distances_mean
    end


    return
end

f, bounds, opt = Metaheuristics.TestProblems.get_problem(:sphere);
optimize(f, bounds, ECA())

```


## Restarting search

```@docs
Restart
```

## Combinatorial Optimization Examples

### Traveling Salesman Problem (TSP)

The TSP seeks the shortest route visiting all cities exactly once and returning to the starting city.

```julia
using Metaheuristics

# Distance matrix for 5 cities
distances = [
    0.0  10.0  15.0  20.0  25.0;
    10.0  0.0  35.0  25.0  30.0;
    15.0 35.0   0.0  30.0  20.0;
    20.0 25.0  30.0   0.0  15.0;
    25.0 30.0  20.0  15.0   0.0
]

# Objective: minimize tour length
function tsp_cost(tour)
    n = length(tour)
    total = sum(distances[tour[i], tour[i+1]] for i in 1:n-1)
    total += distances[tour[end], tour[1]]  # return to start
    return total
end

# Solve using GA with order crossover
n_cities = 5
population_size = 100
ga = GA(
    N = population_size,
    initializer = RandomPermutation(N=population_size),
    crossover = OrderCrossover(), # crossover for permutation encoding
    mutation = SlightMutation()
)

result = optimize(tsp_cost, PermutationSpace(n_cities), ga)
best_tour = minimizer(result)
tour_length = minimum(result)

println("Best tour: ", best_tour)
println("Tour length: ", tour_length)
```

### Knapsack Problem

The knapsack problem is a classic combinatorial optimization problem where the goal is to maximize the value of items that can be carried in a knapsack with a limited capacity. For example, if you have a knapsack with a capacity of 7 kg and you have the following items in the table:

| Item | Profit | Weight (kg) |
|------|--------|--------|
| gold | 1000 | 1 |
| silver | 500 | 2 |
| bronze | 250 | 4 |
| laptop | 800 | 3 |
| camera | 600 | 2 |
| necklace | 300 | 1 |
| watch | 400 | 1 |

The knapsack problem is to find the items that maximize the profit while respecting the capacity constraint.

```julia
using Metaheuristics

# define items and capacity
items = Dict(
    "gold" => (profit=1000, weight=1),
    "silver" => (profit=500, weight=2),
    "bronze" => (profit=250, weight=4),
    "laptop" => (profit=800, weight=3),
    "camera" => (profit=600, weight=2),
    "necklace" => (profit=300, weight=1),
    "watch" => (profit=400, weight=1)
)
capacity = 7

names = collect(keys(items))
profit = [items[k].profit for k in names]
weight = [items[k].weight for k in names]

# load the test problem using the provided parameters
# by default, the problem is encoded as a permutation
f, search_space, optimum = Metaheuristics.TestProblems.knapsack(profit, weight, capacity)

population_size = 100
ga = GA(
    N = population_size,
    initializer = RandomPermutation(N=population_size),
    crossover = OrderCrossover(), # crossover for permutation encoding
    mutation = SlightMutation()
)

# Solve with GA
result = optimize(f, search_space, ga)
selected_items = minimizer(result)[1:findfirst(w -> w > capacity, cumsum(weight[minimizer(result)]))-1]

# table with the selected items
println("Item \t Weight \t Value")
for i in selected_items
    println("$(names[i]) \t $(weight[i]) \t $(profit[i])")
end
println("Total value: ", -minimum(result))  # negative because we minimize


```


## Comparing Multiple Algorithms

Compare performance of different algorithms on the same problem:

```julia
using Metaheuristics

# Define problem
f(x) = sum(x.^2)
bounds = boxconstraints(lb = -10ones(5), ub = 10ones(5))

# Test multiple algorithms
algorithms = [
    ("ECA", ECA()),
    ("DE", DE()),
    ("PSO", PSO()),
    ("ABC", ABC())
]

println("Algorithm Performance Comparison")
println("-" ^ 50)

for (name, algo) in algorithms
    result = optimize(f, bounds, algo)
    println("$name:")
    println("  Minimum: $(minimum(result))")
    println("  F-calls: $(result.f_calls)")
    println("  Time: $(result.final_time - result.start_time) s")
    println()
end
```

## Custom Stopping Criteria Example

Implement diversity-based stopping:

```julia
using Metaheuristics
import LinearAlgebra: norm

# Overwrite stop criteria for ECA
function Metaheuristics.stop_criteria!(
        status,
        parameters::ECA,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    
    if status.stop
        return
    end
    
    # Calculate population diversity
    if length(status.population) > 1
        positions = [sol.x for sol in status.population]
        center = sum(positions) / length(positions)
        avg_distance = sum(p -> norm(p - center), positions) / length(positions)
        
        # Stop when population converges
        if avg_distance < 1e-4
            status.stop = true
            status.stop_msg = "Low population diversity"
        end
    end
    
    return
end

# Now use ECA with custom stopping
f(x) = sum(x.^2)
bounds = boxconstraints(lb = -5ones(10), ub = 5ones(10))
result = optimize(f, bounds, ECA())
println(result.stop_msg)
```
