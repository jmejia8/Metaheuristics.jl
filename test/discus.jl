using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "Discus" begin
    function test_result(result::Vector, fitness::Float64, D::Int, tol::Float64)
        @test ≈(fitness, 0.0, atol=tol)
    end

    # Dimension
    D = 10

    # Objective function
    discus(x::Vector{Float64}) = 1e6x[1].^2 + sum(x[2:end] .^2)

    # ECA results
    result, fitness = eca(discus, D; showResults=false)
    test_result(result, fitness, D, 1e-5)

    # ED results
    result, fitness = DE(discus, D; F = 1, CR = 0.5, showResults=false)
    test_result(result, fitness, D, 1e-5)

    # PSO results
    result, fitness = pso(discus, 2; showResults=false)
    test_result(result, fitness, 2, 1e1)

    # ABC results
    result, fitness = ABC(discus, [-10.0ones(5) 10.0ones(5)]', limit=20)
    test_result(result, fitness, D, 1e-5)

end