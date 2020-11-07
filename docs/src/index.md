# Introduction

Documentation for the Metaheuristics Package.

## Installation

### Julia 0.7 or Later

Open the Julia REPL and press `]` to open the Pkg prompt. To add this package, use the add command:
```
pkg> add https://github.com/jmejia8/Metaheuristics.jl.git
```

## Quick Start

Assume you want to solve the following minimization problem.

Minimize:

$f(x) = \sum_{i=1}^{10} x_i^2$

where $x\in[-10, 10]^{10}$, i.e., $-10 \leq x_i \leq 10$ for $i=1,\ldots,10$.

### Solution

Firstly, import the Metaheuristics package:

```julia
using Metaheuristics
```

Code the objective function:
```julia
f(x) = sum(x.^2)
```

Instantiate the bounds, note that `bounds` should be a $2\times 10$ `Matrix` where
the first row corresponds to the lower bounds whilst the second row corresponds to the
upper bounds.

```julia
bounds = [-10ones(10) 10ones(10)]'
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

