using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "constrained" begin
    function test_result(result::Vector, fitness::Float64, D::Int, tol::Float64)
        @test ≈(fitness, 0.0, atol=tol)
    end

    # Dimension
    D = 5

    # Objective function
    f(x) = (sum((x .- 1).^2), [sum((x .- 1).^2) - 4, sum(sin.(x .- 1)) - 1], [(x[1] - x[2])^2])

    bounds = Array([-10.0ones(D) 10.0ones(D)]')

    # ECA results
    status = optimize(f, bounds, ECA(ε = 1.0, options=Options(debug=false)))

    display(status.best_sol)
    println("--------------------------")
    for ind = status.population
        @show ind.x
    end

    println("==========")
    @show f(ones(D))
    result = status.best_sol.x
    fitness = status.best_sol.f
    test_result(result, fitness, D, 1e-5)

end

