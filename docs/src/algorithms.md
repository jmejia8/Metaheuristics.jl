# Algorithms

List of implemented metaheuristics. The algorithms were implemented based on the
contributor's understanding of the algorithms detailed in the published paper.

| Algorithm | Objective | Constraints | Large Scale | Batch Evaluation| Structure Name         |
|-----------|:-----------|:-----------:|:-----------:| :----------:|------------------------|
| ECA       |  Single    |      ✅     |   ➖        |   ✅        | [`ECA`](@ref)          |
| DE        |  Single    |      ✅     |   ➖        |   ✅        | [`DE`](@ref)           |
| PSO       |  Single    |      ✅     |   ➖        |   ✅        | [`PSO`](@ref)          |
| ABC       |  Single    |      ❌     |   ➖        |   ❌        | [`ABC`](@ref)          |
| MOEA/D-DE |  Multi     |      ➖     |   ➖        |   ❌        | [`MOEAD_DE`](@ref)     |
| GSA       |  Single    |      ❌     |   ❌        |   ✅        | [`CGSA`](@ref)         |
| SA        |  Single    |      ✅     |   ➖        |   ❌        | [`SA`](@ref)           |
| WOA       |  Single    |      ✅     |   ➖        |   ✅        | [`WOA`](@ref)          |
| NSGA-II   |  Multi     |      ✅     |   ➖        |   ✅        | [`NSGA2`](@ref)        |
| NSGA-III  |  Many      |      ✅     |   ➖        |   ✅        | [`NSGA3`](@ref)        |
| SMS-EMOA  |  Multi     |      ✅     |   ➖        |   ✅        | [`SMS_EMOA`](@ref)     |
| SPEA2     |  Multi     |      ✅     |   ➖        |   ✅        | [`SPEA2`](@ref)        |
| BCA       |  Bilevel   |      ✅     |   ❌        |   ❌        | [`BCA`](https://jmejia8.github.io/BilevelHeuristics.jl/dev/algorithms/#BCA) |
| MCCGA     |  Single    |      ❌     |   ❌        |   ❌        | [`MCCGA`](@ref)        |
| GA        |  Single    |      ✅     |   ➖        |   ✅        | [`GA`](@ref)        |
| CCMO      |  Multi     |      ✅     |   ➖        |   ✅        | [`CCMO`](@ref)        |


✅ = supported,
❌ = not supported,
➖ = can be supported by changing default parameters.

- **Batch Evaluation** = Simultaneous evaluation of multiple solutions (batch) see "[Batch Evaluation](@ref)".
- **Constraints** = Equality and inequality constraints.
- **Large Scale** = High dimensional problems (variables space).

## Evolutionary Centers Algorithm

ECA was proposed for solving global optimization problems. See [MejiaMezura2019](@cite) for more information.
```@docs
ECA
```

## Differential Evolution

DE is an evolutionary algorithm based on vector differences.
See [Price2013](@cite) for more details.

```@docs
DE
```

## Particle Swarm Optimization

PSO is a population-based optimization technique inspired by the motion of bird flocks and schooling fish by [KennedyEberhart1995](@cite).

```@docs
PSO
```

## Artificial Bee Colony

A powerful and efficient algorithm for numerical function optimization: artificial bee colony (ABC) algorithm by [KarabogaBasturk2007](@cite).
```@docs
ABC
```


## MOEA/D-DE

Multiobjective optimization problems with complicated Pareto sets by [LiZhang2008](@cite).

```@docs
MOEAD_DE
```


## Gravitational Search Algorithm

Chaotic gravitational constants for the gravitational search algorithm by
[MirjaliliGandomi2017](@cite)

```@docs
CGSA
```


## Simulated Annealing

Physics-inspired algorithm for optimization by [Van1987](@cite).

```@docs
SA
```

## Whale Optimization Algorithm

The Whale Optimization Algorithm inspired by humpback whales proposed in [MirjaliliLewis2016](@cite).

```@docs
WOA
```


## NSGA-II

A fast and elitist multiobjective genetic algorithm: NSGA-II by [Deb2002](@cite).

```@docs
NSGA2
```


## NSGA-III

An Evolutionary Many-Objective Optimization Algorithm Using Reference-Point-Based
Nondominated Sorting Approach, Part I: Solving Problems With Box Constraints by [DebJain2014](@cite).
```@docs
NSGA3
```

## SMS-EMOA

An EMO algorithm using the hypervolume measure as a selection criterion by [Emmerich2005](@cite).
```@docs
SMS_EMOA
```


## SPEA2

Improved strength Pareto evolutionary algorithm by [Zitzler2001](@cite).
```@docs
SPEA2
```

## BCA

Bilevel Centers Algorithm has been proposed to solve bilevel optimization problems.
See [`BilevelHeuristics.BCA`](https://jmejia8.github.io/BilevelHeuristics.jl/dev/algorithms/#BCA) for 
details.


## MCCGA

Machine-coded Compact Genetic Algorithms for real-valued optimization problems by [SatmanAkadal2020mcga](@cite).

```@docs
MCCGA
```

## GA


```@docs
GA
```

## CCMO

A Coevolutionary Framework for Constrained Multiobjective Optimization Problems
proposed by [Tian2020](@cite).

```@docs
CCMO
```

## $\varepsilon$DE

$\varepsilon$ Constrained Differential Evolution with Gradient-Based Mutation and Feasible Elites by [Takahama2006Constrained](@cite).

```@docs
εDE
```
