# Getting Started

After reading this tutorial you'll become an expert in using the `Metaheuristics` module.

## Minimization Problem

Assume you want to optimize the following minimization problem:

Minimize:

$f(x) = 10D + \sum_{i=1}^D x_i^2 - 10cos(2\pi x_i)$

where $x\in [-5, 5]^D$, that is, each coordinate in $x$ is between -5 and 5. Use $D=10$.

Note that the global optimum is obtained when $x_i = 0$ for all $i$. Thus, $\min f(x) = 0$.

Objective function:

```julia
f(x) = 10length(x) + sum( x.^2 - 10cos.(2pi*x) )
```

Bounds:

```julia
bounds = [-5ones(10) 5ones(10)]'
```

## Providing Information

Since the optimum is known, then we can provide this information to the optimizer.

```julia
information = Information(f_optimum = 0.0)
```

## Common Settings

Usually, users could require to limit the number of generations/iterations or the number
of function evaluations. To do that, let's assume that the metaheuristic should evaluate at most
$9000D$ times the objective function. Moreover, since `information` is provided, then we
can set the desired accuracy ($|f(x) - f(x^*)| $) to $10^{-5}$.

```julia
options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)
```

## Choose a Metaheuristic

`Metaheuristics.jl` provides different metaheuristics for optimization such as
Evolutionary Centers Algorithm (ECA), Differential Evolution (DE), Particle Swarm
Optimization (PSO), etc. In this tutorial, we will use `ECA`, but you can use another
algorithm following the same steps.

The metaheuristics accept their parameters but share two common and **optional** settings
`information` and `options`.

```julia
algorithm = ECA(information = information, options = options)
```

!!! warning "Consider changing the default parameters."
    Change population size for high dimensional problems.

## Optimize

Now, we are able to approximate the optimum. To do that it is necessary to use the `optimize`
function as follows:

```julia
result = optimize(f, bounds, algorithm)
```

## Get the Results

Once `optimize` stops, we can get the approximate solutions.

Approximated minimum:

```julia
fx = minimum(result)
```

Approximated minimizer:

```julia
x = minimizer(result)
```

## Get Information about the Resulting Population

Sometimes it is useful to analyze the resulting population (for population-based metaheuristics).
To do that you can use `fvals` to get objective function evaluation and `positions` to
get their positions.

## Bonus

We recommend you wrap your program in a function for performance purposes:

```julia
using Metaheuristics

function main()
    # objective function
    f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )
    
    # limits/bounds
    bounds = [-5ones(10) 5ones(10)]'
    
    # information on the minimization problem
    information = Information(f_optimum = 0.0)

    # generic settings
    options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)
    
    # metaheuristic used to optimize
    algorithm = ECA(information = information, options = options)

    # start the minimization process
    result = optimize(f, bounds, algorithm)

    
    fx = minimum(result)
    x = minimizer(result)

    @show fx
    @show x
end

```

## Summary

Now you are able to approximate global optimum solutions using Metaheuristics.
