using Metaheuristics
using Test



@testset "Box-constrained" begin


    function run_methods(problem)

        desired_accuracy = 1e-4

        f, bounds = Metaheuristics.Benchmark.get_problem(problem)

        information = Information(f_optimum = 0.0)
        options = Options(f_tol = desired_accuracy, seed = 1)

        methods = [
                   ABC(options = options, information = information),
                   CGSA(options = options, information = information),
                   DE(CR = 0.5, options = options, information = information),
                   ECA(options = options, information = information),
                   PSO(options = options, information = information),
                   SA(options = options, information = information),
                   WOA(options = options, information = information),
                  ]

        for method in methods
            fitness = minimum( optimize(f, bounds, method) ) 
            test_result(fitness, desired_accuracy)
        end
    end

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    for problem in [:sphere]
        run_methods(problem)
    end

end
