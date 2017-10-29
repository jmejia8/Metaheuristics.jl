using Metaheuristics
using Base.Test

srand(31415926534)

@testset "Metaheuristics" for tests in [
            "sphere.jl",
            "rastrigin.jl",
            "discus.jl"
]
    include(tests)
end