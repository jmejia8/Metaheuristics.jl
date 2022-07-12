---
title: 'Metaheuristics: A Julia Package for Single- and Multi-Objective Optimization'
tags:
  - Julia
  - optimizaton
  - multi-objective
  - metaheuristics
  - evolutionary algorithms
authors:
  - name: Jesús-Adolfo Mejía-de-Dios
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
  - name: Efrén Mezura-Montes
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Artificial Intelligence Research Institute, University of Veracruz, MEXICO
   index: 1

date: 12 July 2022
bibliography: paper.bib

---

# Summary

`Metaheuristics` is a `Julia` package that implements metaheuristic algorithms for global optimization problems that can contain constraints, and either single or multiple objectives.
Besides, it provides challenging test problems, and performance indicators to assess the performance of the algorithms as usual in Evolutionary Computing.
`Metaheuristics` implements an easy-to-use API that let testing among the implemented algorithm to an interest problem without extra configurations.
The main features and some usage examples are presented in this paper.


# Statement of need



Real-world problems require sophisticated methodologies to provide feasible and efficient solutions of hard problems.
Metaheuristics are algorithms proposed to approximate optimal solutions in a short time to optimization problems unknown mathematical properties. Metaheuristics package implements  state-of-the-art algorithms for constrained, multi-, and many-objective optimization.
Besides, many other utility function have implemented in the package, e.g., performance indicators for performance assessment, scalable benchmark problems, constraint handling techniques, and multi-criteria decision-making methodologies.
To the best knowledge of the authors, `Metaheuristics` is first package in `Julia` containing those ready-to-use features with uniform API.


# Main Features

`Metaheuristics` implements a consistent and intuitive API to perform the approximation of
an objective function bounded in a region which can contains equality or inequality constraints.

$$\min_{x\in X} f(x)$$
subject to
$$
g_j(x)  \leq 0,\ j = 1,2,\ldots,J;
$$
$$
h_l(x)  = 0,\ l = 1,2,\ldots,L;
$$
with $x_{i,\min} \leq x_i \leq x_{i,\max}$, i.e., $X = \prod_{i=1}^D [x_{i,\min},x_{i,\max}]$.

Note that $f$ can be either single- or multi-objective. The main features are listed below:

- A lot of optimization `methods` are implemented.
- Include constraint handling strategies.
- Support for single-, multi-, many-objective optimization. 
- Contain scalable benchmark problems for testing metaheuristics.
- Visualization in console via [UnicodePlots](https://github.com/JuliaPlots/UnicodePlots.jl).
- Performance Indicators for performance assessment.
- Multi-criteria decision-making.


## Example: Single-Objective Optimization

Assume you want to optimize the following minimization problem.
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

It is worth mentioning that `optimize(f, bounds, OPTIMIZER)` where `OPTIMIZER` can be selected from the implemented metaheuristics described in following section.



## Implemented Metaheuristics

List of implemented metaheuristics. The algorithms were implemented based on the
contributor's understanding of the algorithms detailed in the published paper.

Algorithms for single-objective optimization: Evolutionary Centers Algorithms (ECA) [@MejiaMezura2019], Differential evolutionary (DE) [@Price2013], Particle Swarm Optimization (PSO) [@KennedyEberhart1995], Artificial Bee Colony (ABC) [@KarabogaBasturk2007], Gravitational Search Algorithm (GSA) [@MirjaliliGandomi2017], Simulated Annealing (SA) [@Van1987], Whale Optimization Algorithm (WOA) [@MirjaliliLewis2016], generic and compact genetic algorithms [@SatmanAkadal2020mcga].
Besides, multi-objective optimization algorithms such as Multi-objective Evolutionary Algorithm (MOEA/D-DE) [@LiZhang2008], Non-dominated Sorting Genetic Algorithm (NSGA-II and III)[@Deb2002; @DebJain2014], $S$-Metric Selection Evolutionary Multi-objective Algorithm  (SMS-EMOA) [@Emmerich2005], Improved Strength Pareto Evolutionary Algorithm  (SPEA2) [@Zitzler2001],  Coevolutionary Framework for Constrained Multiobjective Optimization (CCMO) [@Tian2020].

## Performance Indicators

Performance Indicators are used to assess algorithms regarding the quality corresponding outcomes as detailed in [@Zitzler2003].
Generational Distance (GD), Inverted Generational Distance (IGD), IGD+, Covering Indicator (C-metric)
Hypervolume (HV), Averaged Hausdorff distance (also known as $\Delta_p$), Spacing Indicator,
and the $\varepsilon$-indicator.

## Multi-Criteria Decision-Making

Multi-Criteria Decision-Making (MCDM) has to be performed after an optimizer reported a set of non-dominated solutions. The implemented MCDM techniques include Compromise Programming [@Ringuest1992] and Region of Interest Archiving [@sebastian2022efficient].
Moreover, an interface for `JMcDM` [@Satman2021JMcDM] which is a package for MCDM that implements many MCDM techniques interfaced within `Metaheuristics`.


# Acknowledgements

The first author acknowledges support from the Mexican Council for Science and Technology (CONACyT) through a scholarship to pursue graduate studies at the University of Veracruz, MEXICO.

# References
