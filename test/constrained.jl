using Metaheuristics
using Test

@testset "Constrained" begin

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    function run_methods_cop()
        desired_accuracy = 1e-4

        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(:constrained1)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        information = Information(f_optimum = 0.0)
        options = Options(f_tol = desired_accuracy, h_tol=1e-5, seed = 2)

        methods = [
                   ECA(options = options, information = information),
                   #DE(options = options, information = information),
                  ]

        for method in methods
            f_calls = 0
            res = optimize(f, bounds, method)
            fitness = minimum( res ) 
            test_result(fitness, desired_accuracy)
            @test f_calls == res.f_calls
        end
    end

    run_methods_cop()
end
