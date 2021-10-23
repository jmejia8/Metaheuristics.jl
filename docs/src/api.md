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
SBX_crossover
```


```@docs
ECA_operator
```

```@docs
DE_crossover
```

```@docs
polynomial_mutation!
```

```@docs
DE_mutation
```

```@docs
MOEAD_DE_reproduction
```

```@docs
binary_tournament
```


```@docs
GA_reproduction
```


```@docs
GA_reproduction_half
```

## Population


```@docs
get_best
```

```@docs
argworst
```

```@docs
argbest
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
non_dominated_sort
```

```@docs
get_fronts
```

```@docs
fast_non_dominated_sort!
```

```@docs
get_non_dominated_solutions_perm
```


```@docs
get_non_dominated_solutions
```


