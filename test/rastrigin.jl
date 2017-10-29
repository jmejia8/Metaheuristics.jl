using Metaheuristics
using Base.Test

# write your own tests here
@testset "Rastrigin" begin
    function test_result(result::Vector, fitness::Float64, D::Int, tol::Float64)
        @test ≈(fitness, 0.0, atol=tol)
    end

    # Dimension
    D = 3

    # Objective function
    rastrigin(x::Vector{Float64}, D=length(x)) = 10D+ sum(x.*x - 10cos.(2π*x))

    result, fitness = eca(rastrigin, D; showResults=false)
    test_result(result, fitness, D, 1e-5)


end