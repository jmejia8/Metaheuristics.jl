using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "Rastrigin" begin
    function test_result(result::Vector, fitness::Float64, D::Int, tol::Float64)
        @test ≈(fitness, 0.0, atol=tol)
    end

    # Dimension
    D = 3

    # Objective function
    rastrigin(x::Vector{Float64}, D=length(x)) = 10D+ sum(x.*x - 10cos.(2π*x))

    result, fitness = eca(rastrigin, D; showResults=false, N = 15D)
    test_result(result, fitness, D, 1e-5)

    # ED results
    result, fitness = DE(rastrigin, D; F = 1, CR = 0.5, showResults=false)
    test_result(result, fitness, D, 1e-5)

    # PSO results
    result, fitness = pso(rastrigin, 2; showResults=false)
    test_result(result, fitness, D, 1e-5)

    # ABC results
    result, fitness = ABC(rastrigin, [-1.0ones(D) 1.0ones(D)]', limit=20)
    test_result(result, fitness, D, 1e-5)

end