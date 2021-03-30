using Metaheuristics
using Test



@testset "Box-constrained" begin


    function run_methods(problem)

        desired_accuracy = 1e-4

        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            global f_calls += 1
            ff(x)
        end

        information = Information(f_optimum = 0.0)
        options = Options(f_tol = desired_accuracy, seed = 1)
        options_2 = Options(f_tol = desired_accuracy, seed = 1, store_convergence=true)

        methods = [
                   ABC(options = options, information = information),
                   CGSA(options = options, information = information),
                   DE(CR = 0.5, options = options_2, information = information),
                   ECA(options = options, information = information),
                   PSO(options = options, information = information),
                   SA(options = options, information = information),
                   WOA(options = options, information = information),
                  ]

        for method in methods
            global f_calls = 0
            res = optimize(f, bounds, method)
            fitness = minimum( res ) 
            test_result(fitness, desired_accuracy)

            # also count the mumber of function evaluations
            @test f_calls == res.f_calls
        end
    end

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    for problem in [:sphere]
        run_methods(problem)
    end

end
