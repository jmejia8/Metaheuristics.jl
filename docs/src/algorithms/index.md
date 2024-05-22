# Index

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


✅ = supported,
❌ = not supported,
➖ = can be supported by changing default parameters.

- **Batch Evaluation** = Simultaneous evaluation of multiple solutions (batch) see "[Batch Evaluation](@ref)".
- **Constraints** = Equality and inequality constraints.
- **Large Scale** = High dimensional problems (variables space).

