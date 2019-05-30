using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

@testset "Metaheuristics" for tests in [
            "sphere.jl",
            "discus.jl",
            "rastrigin.jl"
]
    include(tests)
end