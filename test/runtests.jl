using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
    import Random: seed!
    seed!(31415926534)
end

@testset "Metaheuristics" for tests in [
            "common-methods.jl",
            "box-constrained.jl",
            "constrained.jl",
            "multi-objective.jl",
]
    include(tests)
end
