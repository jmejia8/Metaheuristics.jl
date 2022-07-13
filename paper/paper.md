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

`Metaheuristics` is a `Julia` package that implements metaheuristic algorithms for solving global optimization problems that can contain constraints, and either single or multiple objectives.
The package implements an easy-to-use API that let testing without extra configurations among the implemented optimizers, i.e., to `optimize` an objective function $f(x)$ with $x$ in the corresponding `bounds` using optimizer `ECA`, the user can write `optimize(f, bounds, ECA())` to approximate an optimal solution.
Moreover, `Metaheuristics` provides common features required by the Evolutionary Computing community such as challenging test problems, performance indicators, and other interest utility functions devoted to performance improvement.
This paper presents the main features and some usage examples to illustrate the usage of this package for optimization.


# Statement of need



Real-world problems require sophisticated methodologies to provide feasible and efficient solutions to hard problems.
Metaheuristics are algorithms proposed to approximate optimal solutions in a short time to optimize problems with unknown mathematical properties. Metaheuristics package implements state-of-the-art algorithms for constrained, multi-, and many-objective optimization.
Besides, many other utility functions have been implemented in the package, e.g., performance indicators for performance assessment, scalable benchmark problems, constraint handling techniques, and multi-criteria decision-making methodologies.
To the best knowledge of the authors, `Metaheuristics` is the first package in `Julia` containing those ready-to-use features with a uniform API. The package can be used for both academic and industrial purposes due to it is distributed through a flexible license.



# Main Features

This part aims to describe the most important features included in the package.
First, `Metaheuristics` implements a consistent and intuitive API to perform the approximation to optimal solutions of
an objective function defined as follows.

$$\min_{x\in X} f(x)$$
subject to
$$
g_j(x)  \leq 0,\ j = 1,2,\ldots,J;
$$
$$
h_l(x)  = 0,\ l = 1,2,\ldots,L;
$$
with $x_{i,\min} \leq x_i \leq x_{i,\max}$, i.e., $x \in X = \prod_{i=1}^D [x_{i,\min},x_{i,\max}]$.
Moreover, $f$ can be either single- or multi-objective.

## Implemented Metaheuristics

The list of implemented metaheuristics is detailed in **Table 1**.
Algorithms for single-objective optimization: Evolutionary Centers Algorithms (ECA) [@MejiaMezura2019], Differential evolutionary (DE) [@Price2013], Particle Swarm Optimization (PSO) [@KennedyEberhart1995], Artificial Bee Colony (ABC) [@KarabogaBasturk2007], Gravitational Search Algorithm (GSA) [@MirjaliliGandomi2017], Simulated Annealing (SA) [@Van1987], Whale Optimization Algorithm (WOA) [@MirjaliliLewis2016], generic and compact genetic algorithms [@SatmanAkadal2020mcga].
Besides, `Metaheuristics` includes multi-objective optimization algorithms such as Multi-objective Evolutionary algorithm based on decomposition (MOEA/D-DE) [@LiZhang2008], Non-dominated Sorting Genetic Algorithms (NSGA-II,-III)[@Deb2002; @DebJain2014], $S$-Metric Selection Evolutionary Multi-objective Algorithm  (SMS-EMOA) [@Emmerich2005], Improved Strength Pareto Evolutionary Algorithm  (SPEA2) [@Zitzler2001], and Coevolutionary Framework for Constrained Multiobjective Optimization (CCMO) [@Tian2020].

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

Table: Implemented optimizers so far. Here, "$\checkmark$", "$\times$" and "$-$", mean that
feature in the corresponding column is available, unavailable, or can be enabled by changing default parameters.


## Performance Indicators

Performance Indicators are used to assess algorithms regarding the quality of corresponding outcomes [@Zitzler2003].
`Metaheuristics` implements the following indicators:
Generational Distance (GD), Inverted Generational Distance (IGD), IGD+, Covering Indicator (C-metric),
Hypervolume (HV), Averaged Hausdorff distance (also known as $\Delta_p$), Spacing Indicator,
and the $\varepsilon$-indicator.

## Multi-Criteria Decision-Making

Multi-Criteria Decision-Making (MCDM) has to be performed after an optimizer reported a set of non-dominated solutions. The implemented MCDM techniques include Compromise Programming [@Ringuest1992] and Region of Interest Archiving [@sebastian2022efficient].
Moreover, an interface for `JMcDM` [@Satman2021JMcDM] is a package for MCDM that implements many MCDM techniques interfaced within `Metaheuristics`.



# Installation and Usage

`Metaheuristics` can be installed through the Julia package manager by using the add command in the Pkg prompt:

```julia
pkg> add Metaheuristics
```
Or, equivalently, via the `Pkg` API:
```julia
julia> import Pkg
julia> Pkg.add("Metaheuristics")
```

Now, let us consider the following minimization problem.

Minimize:
$$ f(x) = 10D + \sum_{i=1}^D x_i^2 - 10\cos(2\pi x_i)$$
where the box-constraints (`bounds`) are $x\in [-10, 10]^5$.


We can `optimize` the objective function $f$ with $x$ in corresponding `bounds` as follows:

```julia
julia> f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )

julia> bounds = [-10ones(5) 10ones(5)]

julia> result = optimize(f, bounds)
+=========== RESULT ==========+
  iteration: 543
    minimum: 0
  minimizer: [1.9910474007386918e-10,... , 4.058906169667879e-9]
    f calls: 19005
 total time: 0.1169 s
stop reason: Small difference of objective function values.
+============================+
```

Note that `optimize(f, bounds, OPTIMIZER)` is used to approximate an optimal solution, where `OPTIMIZER` can be selected from the implemented metaheuristics (see **Table 1**).

Finally, It is encouraged to read the  [documentation](https://jmejia8.github.io/Metaheuristics.jl/stable/) for more details, options, and examples about the package.




# Acknowledgements

The first author acknowledges support from the Mexican Council for Science and Technology (CONACyT) through a scholarship to pursue graduate studies at the University of Veracruz, MEXICO.

# References
