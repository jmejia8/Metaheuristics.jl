using Metaheuristics
if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "Sphere" begin
    function test_result(result::Vector, fitness::Float64, D::Int, tol::Float64)
        @test ≈(fitness, 0.0, atol=tol)
    end

    # Dimension
    D = 10

    # Objective function
    sphere(x) = sum(x.*x)

    bounds = Array([-100.0ones(D) 100.0ones(D)]')

    # ECA results
    status = optimize(sphere, bounds)
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-5)

    # ED results
    status = optimize(sphere, bounds, DE(CR=0.5))
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-5)

    # PSO results
    status = optimize(sphere, bounds, PSO())
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-5)

    # ABC results
    result, fitness = ABC(sphere, [-10.0ones(5) 10.0ones(5)]')
    test_result(result, fitness, D, 1e-5)
end
