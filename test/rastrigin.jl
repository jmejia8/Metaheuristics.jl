using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
    import Random: seed!
    seed!(31415926534)
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

    bounds = Array([-100.0ones(D) 100.0ones(D)]')


    status = optimize(rastrigin, bounds, ECA())
    result = status.best_sol.x
    fitness = status.best_sol.f

    test_result(result, fitness, D, 1e-5)

    # ED results
    status = optimize(rastrigin, bounds, DE(CR=0.5))
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-5)

    # PSO results
    pso = PSO(N=300)
    pso.options.f_calls_limit = 30000D
    status = optimize(rastrigin, bounds, pso)
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-3)

    # ABC results
    result, fitness = ABC(rastrigin, [-1.0ones(D) 1.0ones(D)]', limit=20)
    test_result(result, fitness, D, 1e-5)

end
