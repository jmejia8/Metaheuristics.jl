using Metaheuristics
using Test

@testset "Constrained" begin

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    desired_accuracy = 1e-4

    f, bounds, optimums = Metaheuristics.TestProblems.get_problem(:constrained1)

    information = Information(f_optimum = 0.0)
    options = Options(f_tol = desired_accuracy, h_tol=1e-5, seed = 1)

    methods = [
               CECA(options = options, information = information),
               # DE(options = options, information = information),
              ]

    for method in methods
        res = optimize(f, bounds, method)
        display(res)
        fitness = minimum( res ) 
        test_result(fitness, desired_accuracy)
    end
end
