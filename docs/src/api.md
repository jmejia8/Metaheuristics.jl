# API References


```@docs
optimize
```

```@docs
State
```

```@docs
Information
```


```@docs
Options
```


```@docs
convergence
```


```@docs
minimizer
```

```@docs
minimum(state::State)
```

```@docs
Metaheuristics.TerminationStatusCode
```


```@docs
termination_status_message
```

## Methods for Solutions/Individuals

```@docs
positions
```

```@docs
fvals
```

```@docs
nfes
```



```@docs
Metaheuristics.create_child
```

```@docs
get_position
```


```@docs
fval
```

```@docs
Metaheuristics.gval
```

```@docs
Metaheuristics.hval
```

```@docs
Metaheuristics.is_feasible
```

```@docs
Metaheuristics.is_better
```


```@docs
Metaheuristics.dominates
```


```@docs
Metaheuristics.compare
```


```@docs
Metaheuristics.gen_initial_state
```

## Variation 

```@docs
Metaheuristics.SBX
```

```@docs
Metaheuristics.SBX_crossover
```


```@docs
Metaheuristics.ECA_operator
```

```@docs
Metaheuristics.DE_crossover
```

```@docs
Metaheuristics.polynomial_mutation!
```

```@docs
Metaheuristics.DE_mutation
```

```@docs
Metaheuristics.MOEAD_DE_reproduction
```

```@docs
Metaheuristics.binary_tournament
```


```@docs
Metaheuristics.GA_reproduction
```


```@docs
Metaheuristics.GA_reproduction_half
```

```@docs
RandomBinary
```

```@docs
RandomPermutation
```

```@docs
TournamentSelection
```
```@docs
RouletteWheelSelection
```

```@docs
UniformCrossover
```
```@docs
OrderCrossover
```

```@docs
BitFlipMutation
```

```@docs
SlightMutation
```


```@docs
PolynomialMutation
```

```@docs
GenerationalReplacement
```

```@docs
ElitistReplacement
```

## Population


```@docs
Metaheuristics.get_best
```

```@docs
Metaheuristics.argworst
```

```@docs
Metaheuristics.argbest
```

```@docs
nadir
```

```@docs
ideal
```


```@docs
Metaheuristics.non_dominated_sort
```

```@docs
Metaheuristics.get_fronts
```

```@docs
Metaheuristics.fast_non_dominated_sort!
```

```@docs
Metaheuristics.get_non_dominated_solutions_perm
```


```@docs
Metaheuristics.get_non_dominated_solutions
```



## Stopping Criteria


```@docs
Metaheuristics.diff_check
```


```@docs
Metaheuristics.call_limit_stop_check
```

```@docs
Metaheuristics.iteration_stop_check
```

```@docs
Metaheuristics.time_stop_check
```

```@docs
Metaheuristics.accuracy_stop_check
```

```@docs
Metaheuristics.var_stop_check
```

## Sampling

```@docs
Metaheuristics.sample
```

```@docs
Metaheuristics.RandomInBounds
```


```@docs
Metaheuristics.LatinHypercubeSampling
```


```@docs
Metaheuristics.Grid
```

