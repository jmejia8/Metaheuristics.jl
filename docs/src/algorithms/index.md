# Index

List of implemented metaheuristics. The algorithms were implemented based on the
contributor's understanding of the algorithms detailed in the published paper.

## Quick Selection Guide

### By Problem Type

- **Box-constrained (continuous)**: Use [`ECA`](@ref), [`DE`](@ref), [`PSO`](@ref), or [`SHADE`](@ref)
- **Constrained optimization**: Use [`εDE`](@ref), [`ECA`](@ref), or [`SHADE`](@ref)
- **Multi-objective (2-3 objectives)**: Use [`NSGA2`](@ref), [`SPEA2`](@ref), or [`SMS_EMOA`](@ref)
- **Many-objective (>3 objectives)**: Use [`NSGA3`](@ref)
- **Constrained multi-objective**: Use [`NSGA2`](@ref), [`NSGA3`](@ref), or [`CCMO`](@ref)
- **Permutation-based**: Use [`GA`](@ref) with [`OrderCrossover`](@ref) or [`BRKGA`](@ref)
- **Binary**: Use [`GA`](@ref) with [`BitFlipMutation`](@ref)
- **Large-scale**: Use [`CSO`](@ref)

### Performance Characteristics

- **Fast convergence**: [`ECA`](@ref), [`SHADE`](@ref)
- **Good for multimodal**: [`DE`](@ref), [`PSO`](@ref), [`ABC`](@ref)
- **Robust**: [`NSGA2`](@ref), [`ECA`](@ref)
- **Simple to use**: [`ECA`](@ref), [`DE`](@ref), [`PSO`](@ref)

### Batch Evaluation

Algorithms marked with ✅ in "Batch Evaluation" column support simultaneous evaluation
of multiple solutions, which is useful for:
- Parallel computing with threads or distributed systems
- GPU acceleration
- Expensive function evaluations

Set `parallel_evaluation=true` in [`Options`](@ref) to enable batch evaluation.
See [Parallelization](@ref) tutorial for details.

## Algorithm Table

| Algorithm | Objective | Constraints | Large Scale | Batch Evaluation| Structure Name         |
|-----------|:-----------|:-----------:|:-----------:| :----------:|------------------------|
| ECA       |  Single    |      ✅     |   ➖        |   ✅        | [`ECA`](@ref)          |
| DE        |  Single    |      ✅     |   ➖        |   ✅        | [`DE`](@ref)           |
| PSO       |  Single    |      ✅     |   ➖        |   ✅        | [`PSO`](@ref)          |
| ABC       |  Single    |      ❌     |   ➖        |   ❌        | [`ABC`](@ref)          |
| MOEA/D-DE |  Multi     |      ➖     |   ➖        |   ❌        | [`MOEAD_DE`](@ref)     |
| GSA       |  Single    |      ❌     |   ❌        |   ✅        | [`CGSA`](@ref)         |
| SA        |  Single    |      ✅     |   ➖        |   ❌        | [`SA`](@ref)           |
| NSGA-II   |  Multi     |      ✅     |   ➖        |   ✅        | [`NSGA2`](@ref)        |
| NSGA-III  |  Many      |      ✅     |   ➖        |   ✅        | [`NSGA3`](@ref)        |
| SMS-EMOA  |  Multi     |      ✅     |   ➖        |   ✅        | [`SMS_EMOA`](@ref)     |
| SPEA2     |  Multi     |      ✅     |   ➖        |   ✅        | [`SPEA2`](@ref)        |
| BCA       |  Bilevel   |      ✅     |   ❌        |   ❌        | [`BCA`](https://jmejia8.github.io/BilevelHeuristics.jl/dev/algorithms/#BCA) |
| MCCGA     |  Single    |      ❌     |   ❌        |   ❌        | [`MCCGA`](@ref)        |
| GA        |  Single    |      ✅     |   ➖        |   ✅        | [`GA`](@ref)           |
| CCMO      |  Multi     |      ✅     |   ➖        |   ✅        | [`CCMO`](@ref)         |
| $\varepsilon$DE |  Single   | ✅     |   ➖        |   ✅        | [`εDE`](@ref)          |
| BRKGA     |  Single    |      ✅     |   ➖        |   ✅        | [`BRKGA`](@ref)           |
| SHADE     |  Single    |      ✅     |   ➖        |   ✅        | [`SHADE`](@ref)           |
| CSO       |  Single    |      ✅     |   ✅        |   ✅        | [`CSO`](@ref)           |


✅ = supported,
❌ = not supported,
➖ = can be supported by changing default parameters.

- **Batch Evaluation** = Simultaneous evaluation of multiple solutions (batch) see "[Batch Evaluation](@ref)".
- **Constraints** = Equality and inequality constraints.
- **Large Scale** = High dimensional problems (variables space).

