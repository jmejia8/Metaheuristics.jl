# Introduction

Optimization is one of the most common task in the scientific and industry field but
real-world problems require high-performance algorithms to optimize non-differentiable,
non-convex, dicontinuous functions. Different metaheuristics algorithms have been
proposed to solve optimization problems but without strong assumptions about the objective
function.

This package implements state-of-the-art metaheuristics algorithms for global optimization.
The aim of this package is to provide easy to use (and fast) metaheuristics for numerical
global optimization.

## Installation

Open the Julia (Julia 0.7 or Later) REPL and press `]` to open the Pkg prompt. To add this package, use the add command:

```
pkg> add https://github.com/jmejia8/Metaheuristics.jl.git
```

## Quick Start

Assume you want to solve the following minimization problem.

![Rastrigin Surface](figs/rastrigin.png)

Minimize:

$f(x) = 10D + \sum_{i=1}^{D}  x_i^2 - 10\cos(2\pi x_i)$

where $x\in[-5, 5]^{D}$, i.e., $-5 \leq x_i \leq 5$ for $i=1,\ldots,D$. $D$ is the
dimension number, assume $D=10$.

### Solution

Firstly, import the Metaheuristics package:

```julia
using Metaheuristics
```

Code the objective function:
```julia
f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x)  )
```

Instantiate the bounds, note that `bounds` should be a $2\times 10$ `Matrix` where
the first row corresponds to the lower bounds whilst the second row corresponds to the
upper bounds.

```julia
D = 10
bounds = [-5ones(D) 5ones(D)]'
```

Approximate the optimum using the function `optimize`.

```julia
result = optimize(f, bounds)
```

Optimize returns a `State` datatype which contains some information about the approximation.
For instance, you may use mainly two functions to obtain such approximation.

```julia
@show minimum(result)
@show minimizer(result)
```
