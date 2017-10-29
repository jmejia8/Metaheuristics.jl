# Metaheuristics

A Julia package for metaheuristic optimization algorithms. Evolutionary are considered.

## Installation

```julia
Pkg.clone("git@github.com:jmejia8/Metaheuristics.jl.git")
```

## Algorithms

- ECA algorithm
- Differential Evolution (DE) algorithm

### ECA

ECA is a new metaheuristic optimization algorithm based on center of mass. ECA minimizes a non-negative objective function.

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
- **limits:** bounds for variables.

#### Example
```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

# Dimension
D = 10

result, fitness = eca(sphere, D; limits=(-10, 10))

```

### DE
Differential Evolution is a method that optimizes a problem by iteratively trying to improve a candidate solution with regard to a given measure of quality. [Read more...](https://en.wikipedia.org/wiki/Differential_evolution)

#### Parameters
- **func:** objective function 
- **D:** dimension.
- **F:** DE-stepsize F_weight from interval [0, 2].
- **N:** Number of population members.
-**CR:** Crossover probability constant from interval [0, 1].
- **max_evals:** number evaluations
- **strategy:** DE strategy
	- `:rand1` DE/rand/1
	- `:rand2` DE/rand/2             
	- `:best1` DE/best/1             
	- `:best2` DE/best/2             
	- `:randToBest1` DE/rand-to-best/1             
- **termination:** criteria function for algorithm termination
- **showResults:** show details of fitness population values
- **limits:** bounds for variables.

#### Example

```julia
using Metaheuristics

# Objective function
sphere(x) = sum(x.^2)

# Dimension
D = 10

result, fitness = diffEvolution(sphere, D; limits=(-10, 10))

```