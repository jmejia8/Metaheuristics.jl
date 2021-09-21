# Examples

This package provides different tools for optimization. Hence, this section gives different
examples for using the implemented `Metaheuristics`.

## Single-Objective Optimization


```@repl
using Metaheuristics

f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) ) # objective function

bounds = [-5ones(10) 5ones(10)]' # limits/bounds

information = Information(f_optimum = 0.0); # information on the minimization problem

options = Options(f_calls_limit = 9000*10, f_tol = 1e-5); # generic settings

algorithm = ECA(information = information, options = options) # metaheuristic used to optimize

result = optimize(f, bounds, algorithm) # start the minimization proccess


minimum(result)
minimizer(result)


result = optimize(f, bounds, algorithm) # note that second run is faster

```

## Providing Initial Solutions

Sometimes you may need to use the starter solutions you need before the optimization
process begins, well, this example illustrates how to do it.

```@repl
using Metaheuristics
f, bounds, optimums = Metaheuristics.TestProblems.get_problem(:sphere);
D = size(bounds,2);

x_known = 0.6ones(D) # known solution

X = [ bounds[1,:] + rand(D).* ( bounds[2,:] -  bounds[1,:]) for i in 1:19  ]; # random solutions (uniform distribution)

push!(X, x_known); # save an interest solution

population = [ Metaheuristics.create_child(x, f(x)) for x in X ]; # generate the population with 19+1 solutions

prev_status = State(Metaheuristics.get_best(population), population); # prior state

method = ECA(N = length(population))
method.status = prev_status; # say to ECA that you have generated a population

optimize(f, bounds, method) # optimize
```


## Constrained Optimization

It is common that optimization models include constraints that must be satisfied for example:
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

!!! compat "Constriants handling"
    In this package, if the algorithm was not designed for constrained optimization,
    then solutions with lower constraint violation sum will be preferred.

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

bounds = [-2.0 -2; 2 2]

optimize(f, bounds, ECA(N=30, K=3))
```

## Multiobjective Optimization

To implement a multiobjective optimization and solve it, you can proceed as usual. Here,
you need to provide constraints if they exist, otherwise put `gx = [0.0]; hx = [0.0];`
to indicate an unconstrained multiobjective problem

```@repl
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

bounds = [zeros(30) ones(30)]';

optimize(f, bounds, NSGA2())
```

## Modifying an Existing Metaheuristic

You may need to modify one of the implemented metaheuristic to improve the
algorithm performance or test new mechanisms. This example illustrate how to do it.


!!! warning "Modifying algorithms could break stuff"
    It is recommended to put the new methods in `module`s rather than in global scope in
    order to avoid errors.


Let's assume that we want to modify the stop criteria for `ECA`. See [Contributing](@ref) 
for more details.

```@repl
using Metaheuristics
import LinearAlgebra: norm

# overwrite method
function Metaheuristics.stop_criteria!(
        status,
        parameters::ECA, # It is important indicate the modified Metaheuristic 
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

    # (optional and not recommended) print when this critaria is met
    if status.stop
        @info "Diversity-based stop criteria"
        @show distances_mean
    end


    return
end

f, bounds, opt = Metaheuristics.TestProblems.get_problem(:sphere);
optimize(f, bounds, ECA())

```
