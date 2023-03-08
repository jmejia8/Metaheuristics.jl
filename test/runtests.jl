using Metaheuristics
using Test
using Aqua # Aqua: Auto QUality Assurance for Julia packages
using UnicodePlots

Aqua.test_all(Metaheuristics, ambiguities = false)

@testset "Ambiguities" begin
    Aqua.test_ambiguities(Metaheuristics, recursive = false)
end

for tests in [
              "common-methods.jl",
              "box-constrained.jl",
              "optimize_api.jl",
              "constrained.jl",
              "multi-objective.jl",
              "combinatorial.jl",
              "decisionmaking.jl",
              "initial_solution.jl",
             ]
    include(tests)
end
