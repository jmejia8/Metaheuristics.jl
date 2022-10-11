using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
    import Random: seed!
    seed!(31415926534)
end

for tests in [
              "common-methods.jl",
              "box-constrained.jl",
              "optimize_api.jl",
              "constrained.jl",
              "multi-objective.jl",
              "combinatorial.jl",
              "decisionmaking.jl",
             ]
    include(tests)
end
