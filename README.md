# Metaheuristics

A Julia package for metaheuristic optimization algorithms. Evolutionary are considered.

## Installation

```julia
Pkg.clone("git@github.com:jmejia8/Metaheuristics.jl.git")
```

## Algorithms

- ECA algorithm

### ECA

ECA is a new metaheuristic optimization algorithm based on center of mass. ECA minimizes an non-negative objective function.

#### Parameters
- **mfunc:** objective function 
- **D:** dimension
- **η_max:** stepsize
- **K:** number of neighbors for generating the center of mass.
- **N:** population size
- **max_evals:** number evaluations
- **termination:** criteria function for algorithm termination
- **showResults:** show details of fitness population values
- **correctSol:** if true, it corrects the solution
- **limits:** bound for variables.

#### Example
```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

# Dimension
D = 10

result, fitness = eca(sphere, D; limits=(-10, 10))

```