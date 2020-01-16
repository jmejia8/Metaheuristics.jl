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
            "sphere.jl",
            "discus.jl",
            "rastrigin.jl",
            "constrained.jl"
]
    include(tests)
end