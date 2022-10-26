---
title: 'Metaheuristics: A Julia Package for Single- and Multi-Objective Optimization'
tags:
  - Julia
  - optimization
  - multi-objective
  - metaheuristics
  - evolutionary algorithms
authors:
  - name: Jesús-Adolfo Mejía-de-Dios
    orcid: 0000-0002-0367-2967
    corresponding: true
    affiliation: 1
  - name: Efrén Mezura-Montes
    orcid: 0000-0002-1565-5267 
    affiliation: 1
affiliations:
 - name: Artificial Intelligence Research Institute, University of Veracruz, MEXICO
   index: 1

date: 12 July 2022
bibliography: paper.bib

---

# Summary

Metaheuristics is a Julia package that implements metaheuristic algorithms for solving global optimization problems that can contain constraints and either single or multiple objectives.
The package exposes an easy-to-use API to enable testing without requiring lengthy configuration for choosing among the implemented optimizers. For example, to optimize an objective function $f(x)$ where solution $x$ is within the corresponding bounds using the Evolutionary Centers Algorithm optimizer, the user can call `optimize(f, bounds, ECA())`.
Moreover, Metaheuristics provides common features required by the evolutionary computing community such as challenging test problems, performance indicators, and other notable utility functions.
This paper presents the main features followed by some examples to illustrate the usage of this package for optimization.


# Statement of need

Real-world problems often require sophisticated methods to solve.
Metaheuristics are stochastic algorithms that approximate optimal solutions quickly where not all of the mathematical properties of the problem are known.
Our package implements state-of-the-art algorithms for constrained, multi-, and many-objective optimization.
It also includes many other utility functions such as performance indicators, scalable benchmark problems, constraint handling techniques, and multi-criteria decision-making methodologies.
Although similar software has been proposed in different programming languages such as Python [@pymoo], MATLAB [@PlatEMO], C/C++ [@Biscani2020], and Java [@jmetal], among others [@nlopt],
Metaheuristics is the first package in Julia containing ready-to-use metaheuristic algorithms and utility functions with a uniform API.

There are also some packages implemented in Julia for global optimization. For example, `Optim.jl` [@mogensen2018optim] implements global optimizers such as Particle Swarm Optimization, `BlackBoxOptim.jl` [@RobertFeldtBBO] implements a couple of stochastic heuristics for black-box optimization, `Evolutionary.jl` [@artewilde2021] is a framework for evolutionary computing, and `CMAEvolutionStrategy.jl` [@CMAEvolutionStrategy] implements a CMA Evolution Strategy. Unlike these packages, Metaheuristics can handle equality and inequality constraints and supports multi-objective problems.

# Main Features

This section describes the primary features of Metaheuristics.
First, Metaheuristics implements a consistent and intuitive API to approximate the optimal solutions of the following minimization problem:

$$\min_{x\in X} f(x)$$
subject to
$$
g_j(x)  \leq 0,\ j = 1,2,\ldots,J;
$$
$$
h_l(x)  = 0,\ l = 1,2,\ldots,L;
$$
where $x_{i,\min} \leq x_i \leq x_{i,\max}$, i.e., $x \in X = \prod_{i=1}^D [x_{i,\min},x_{i,\max}]$.
$f$ can be either single- or multi-objective.

## Implemented Metaheuristics

Implemented metaheuristic algorithms are detailed in **Table 1**.
Algorithms for single-objective optimization include Evolutionary Centers Algorithm (ECA) [@MejiaMezura2019], Differential Evolution (DE) [@Price2013], Particle Swarm Optimization (PSO) [@KennedyEberhart1995], Artificial Bee Colony (ABC) [@KarabogaBasturk2007], Gravitational Search Algorithm (GSA) [@MirjaliliGandomi2017], Simulated Annealing (SA) [@Van1987], Whale Optimization Algorithm (WOA) [@MirjaliliLewis2016], and Machine Coded Compact Genetic Algorithms (MCCGA) [@SatmanAkadal2020mcga].
Metaheuristics also includes multi-objective optimization algorithms such as a Multi-Objective Evolutionary Algorithm Based on Decomposition (MOEA/D-DE) [@LiZhang2008], Non-dominated Sorting Genetic Algorithms (NSGA-II,-III)[@Deb2002; @DebJain2014], $S$-Metric Selection Evolutionary Multi-objective Algorithm  (SMS-EMOA) [@Emmerich2005], Improved Strength Pareto Evolutionary Algorithm (SPEA2) [@Zitzler2001]m and Coevolutionary Framework for Constrained Multiobjective Optimization (CCMO) [@Tian2020].

| Algorithm | Objective  | Constraint Handling |  Batch Evaluation     | Authors         |
|---------------|:--------|:----------:|:------------:|:---------------------------|
| ECA       |  Single | $\checkmark$ | $\checkmark$ |  @MejiaMezura2019  |
| DE        |  Single | $\checkmark$ | $\checkmark$ |  @Price2013  |
| PSO       |  Single | $\checkmark$ | $\checkmark$ |  @KennedyEberhart1995  |
| ABC       |  Single | $\times$     | $\times$     |  @KarabogaBasturk2007  |
| GSA       |  Single | $\times$     | $\checkmark$ |  @MirjaliliGandomi2017   |
| SA        |  Single | $\checkmark$ | $\times$     |  @Van1987   |
| WOA       |  Single | $\checkmark$ | $\checkmark$ |  @MirjaliliGandomi2017   |
| MCCGA     |  Single | $\times$     | $\times$     |  @SatmanAkadal2020mcga   |
| GA        |  Single | $\checkmark$ | $\checkmark$ |  @goldberg2002design   |
| MOEA/D-DE |  Multi  | $-$          | $\times$     |  @LiZhang2008   |
| NSGA-II   |  Multi  | $\checkmark$ | $\checkmark$ |  @Deb2002   |
| SMS-EMOA  |  Multi  | $\checkmark$ | $\checkmark$ |  @Emmerich2005   |
| SPEA2     |  Multi  | $\checkmark$ | $\checkmark$ |  @Zitzler2001   |
| CCMO      |  Multi  | $\checkmark$ | $\checkmark$ |  @Tian2020   |
| NSGA-III  |  Many   | $\checkmark$ | $\checkmark$ |  @DebJain2014   |

Table: Implemented optimizers so far. Here, "$\checkmark$", "$\times$" and "$-$", respectively mean that
the feature in the corresponding column is available, unavailable, or can be enabled by changing default parameters.

## Performance Indicators

Performance indicators are used to assess algorithms with respect to the quality of the corresponding outcomes [@Zitzler2003].
Metaheuristics implements Generational Distance (GD), Inverted Generational Distance (IGD), IGD+, Covering Indicator (C-metric),
Hypervolume (HV), Averaged Hausdorff distance (also known as $\Delta_p$), Spacing Indicator, and the $\varepsilon$-indicator.

## Multi-Criteria Decision-Making

Metaheuristics also provides algorithms for when Multi-Criteria Decision-Making (MCDM) has to be performed after an optimizer has reported a set of non-dominated solutions. The implemented MCDM techniques include Compromise Programming [@Ringuest1992] and Region of Interest Archiving [@sebastian2022efficient].
Metaheuristics also includes an interface to the JMcDM package [@Satman2021JMcDM], which implements many further MCDM techniques.

# Installation and Usage

Metaheuristics can be installed through the Julia package manager by using the add command in the Pkg prompt:

```julia
pkg> add Metaheuristics
```
Or, equivalently, via the `Pkg` API:
```julia
julia> import Pkg
julia> Pkg.add("Metaheuristics")
```

Now, let us consider the following minimization problem:

Minimize:
$$ f(x) = 10D + \sum_{i=1}^D x_i^2 - 10\cos(2\pi x_i)$$
with box constraints (`bounds`) $x\in [-10, 10]^5$.

We can optimize the objective function $f$ with solution $x$ within the desired bounds as follows:

```julia
julia> f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )

julia> bounds = [-10ones(5) 10ones(5)]

julia> result = optimize(f, bounds, DE(CR=0.5))
+=========== RESULT ==========+
  iteration: 550
    minimum: 0
  minimizer: [-1.5469066028117595e-9,..., 3.797322900567224e-9]
    f calls: 27500
 total time: 0.0388 s
stop reason: Small difference of objective function values.
+============================+
```

`optimize(f, bounds, OPTIMIZER)` is used to approximate an optimal solution, where `OPTIMIZER` can be selected from the implemented metaheuristics (see **Table 1**). Note that `OPTIMIZER` is a stochastic procedure and therefore each run may show different outputs. 

Finally, potential users are encouraged to read the [documentation](https://jmejia8.github.io/Metaheuristics.jl/stable/) for more details, options, and examples.

# Acknowledgements

The first author acknowledges support from the Mexican Council for Science and Technology (CONACyT) through a scholarship to pursue graduate studies at the University of Veracruz, MEXICO.

# References
