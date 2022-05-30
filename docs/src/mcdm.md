# Multi-Criteria Decision Making 

A set of Multi-Criteria Decision Making (MCDM) methods are available in `Metaheuristics.jl`.

!!! compat "Maximization or Minimization"
    Here, minimization is always assumed.   

Firstly, it is recommended to read the following two functions.

```docs
decisionmaking
```


```docs
best_alternative
```

## JMcDM


[JMcDM](https://github.com/jbytecode/JMcDM) is a package for MCDM developed by [Satman2021JMcDM](@cite). Many methods have been implemented there, and many of them have been interfaced here.

The Main method to use JMcDM in Metaheuristics is described as follows.

```julia
mcdm(fs, w, method)
```

Perform selected `method` for a given `fs` and weight vector `w`.
Here, `fs` can be a set of non-dominated solutions (population), a `State`
or a decision matrix.

Also, `method` can be selected from [JMcDM](https://github.com/jbytecode/JMcDM) package.

Supported MCDM methods:

* `ArasMethod`
* `CocosoMethod`
* `CodasMethod`
* `CoprasMethod`
* `EdasMethod`
* `ElectreMethod`
* `GreyMethod`
* `MabacMethod`
* `MaircaMethod`
* `MooraMethod`
* `SawMethod`
* `TopsisMethod`
* `VikorMethod`
* `WPMMethod`
* `WaspasMethod`
* `MarcosMethod`

See the [JMcDM documentation](https://jbytecode.github.io/JMcDM/docs/build/mcdms/) for
more details about the methods.

### Example 1:

Performing MCDM using a population.

```julia-repl
julia> using JMcDM

julia> _, _, population = Metaheuristics.TestProblems.ZDT1();

julia> dm = mcdm(population, [0.5, 0.5], TopsisMethod());

julia> population[dm.bestIndex]
(f = [0.5353535353535354, 0.2683214262030523], g = [0.0], h = [0.0], x = [5.354e-01, 0.000e+00, …, 0.000e+00])
```

### Example 2:

Performing MCDM using results from metaheuristic.

```julia-repl
julia> using JMcDM

julia> f, bounds, _ = Metaheuristics.TestProblems.ZDT1();

julia> res = optimize(f, bounds, NSGA2());

julia> dm = mcdm(res, [0.5, 0.5], TopsisMethod());

julia> res.population[dm.bestIndex]
(f = [0.32301132058506055, 0.43208538139854685], g = [0.0], h = [0.0], x = [3.230e-01, 1.919e-04, …, 1.353e-04])
```

### Selecting best alternative

```julia
best_alternative(res, w, method)
```
Perform McDM using results from metaheuristic and return best alternative in `res.population`.


```julia-repl
julia> f, bounds, _ = Metaheuristics.TestProblems.ZDT1();

julia> res = optimize(f, bounds, NSGA2());

julia> best_sol = best_alternative(res, [0.5, 0.5], TopsisMethod())
(f = [0.32301132058506055, 0.43208538139854685], g = [0.0], h = [0.0], x = [3.230e-01, 1.919e-04, …, 1.353e-04])
```


##  Region of Interest Archiving

```@docs
ROIArchiving
```

## Compromise Programming

```docs
CompromiseProgramming
```
