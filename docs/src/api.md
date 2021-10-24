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
pareto_front(st::State)
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
gval
```

```@docs
hval
```

```@docs
is_feasible
```

```@docs
is_better
```


```@docs
dominates
```


```@docs
compare
```

## Variation Operators

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
pareto_front(st::Array)
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


