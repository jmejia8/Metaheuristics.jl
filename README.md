# Metaheuristics

A Julia package for metaheuristic optimization algorithms. Evolutionary are considered.

[![Build Status](https://travis-ci.org/jmejia8/Metaheuristics.jl.svg?branch=master)](https://travis-ci.org/jmejia8/Metaheuristics.jl)
[![Coverage Status](https://coveralls.io/repos/github/jmejia8/Metaheuristics.jl/badge.svg?branch=master)](https://coveralls.io/github/jmejia8/Metaheuristics.jl?branch=master)

## Installation

### Julia 0.7 or Later

Open the Julia REPL and press `]` to open the Pkg prompt. To add this package, use the add command:
```
pkg> add https://github.com/jmejia8/Metaheuristics.jl.git
```

## Algorithms

- ECA algorithm
- Differential Evolution (DE) algorithm
- Particle swarm optimization (PSO) algorithm

### Optimize

`optimize` function is used to optimize a D-dimensional function: `optimize(f::Function, bounds::Array, method::AbstractAlgorithm )`

- `f` objective function
- `bounds` a 2 times D matrix that contains the lower and upper bounds by rows.
- `method` optimization method: `ECA`, `DE`.

### ECA

ECA() is a new metaheuristic optimization algorithm based on center of mass. ECA minimizes an objective function.

#### Parameters
- `η_max:` stepsize.
- `K:` number of neighbors for generating the center of mass.
- `N:` population size.

#### Example
```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

bounds = [-10 -10 -10 -10;
             10  10  10  10
]

eca = ECA()

result = optimize(sphere, bounds, eca)

```

### DE
Differential Evolution `DE` is a method that optimizes a problem by iteratively trying to improve a candidate solution with regard to a given measure of quality. [Read more...](https://en.wikipedia.org/wiki/Differential_evolution)

#### Parameters
- `F:` DE-stepsize F_weight from interval [0, 2].
- `N:` Number of population members.
- `CR:` Crossover probability constant from interval [0, 1].
- `strategy:` DE strategy
	- `:rand1` DE/rand/1
	- `:rand2` DE/rand/2             
	- `:best1` DE/best/1             
	- `:best2` DE/best/2             
	- `:randToBest1` DE/rand-to-best/1             

#### Example

```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

# bounds = [-10 -10 -10 -10;
             10  10  10  10
]

de = DE()

result = optimize(sphere, bounds, de)

```


<!-- ### PSO
Particle swarm optimization is a population based stochastic optimization technique developed by Dr. Eberhart and Dr. Kennedy  in 1995, inspired by social behavior of bird flocking or fish schooling. [Read more...](https://en.wikipedia.org/wiki/Particle_swarm_optimization)

#### Parameters
- **func:** objective function.
- **D:** dimension.
- **N:** Number of population members.
- **C1, C2**  learning factors (C1 = C2 = 2).
- **ω:** Inertia weight used for balancing the global search.
- **max_evals:** number evaluations.
- **termination:** criteria function for algorithm termination.
- **showResults:** show details of fitness population values.
- **limits:** bounds for variables.

#### Example

```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

# Dimension
D = 10

result, fitness = pso(sphere, D; limits=(-10, 10))

``` -->