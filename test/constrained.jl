using Metaheuristics
using Test

@testset "Box-constrained" begin

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    desired_accuracy = 1e-4

    f, bounds = Metaheuristics.Benchmark.get_problem(:constrained1)

    information = Information(f_optimum = 0.0)
    options = Options(f_tol = desired_accuracy, seed = 1)

    methods = [
               ECA(options = options, information = information),
               DE(options = options, information = information),
              ]

    for method in methods
        fitness = minimum( optimize(f, bounds, method) ) 
        test_result(fitness, desired_accuracy)
    end
end
