using Metaheuristics
using Test

@testset "Constrained" begin

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    function run_methods_cop(problem)
        desired_accuracy = 1e-4

        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        information = Information(f_optimum = 0.0)
        options = Options(f_tol = desired_accuracy, h_tol=1e-5, seed = 2, debug=false, store_convergence=true)

        methods = [
                   ECA(ε=0.5, options = options, information = information),
                   εDE(options = options, information = information),
                  ]

        for method in methods
            f_calls = 0
            res = optimize(f, bounds, method)
            show(IOBuffer(), res)
            show(IOBuffer(), "text/plain", res.population)
            show(IOBuffer(), "text/html", res.population)
            show(IOBuffer(), res.population[1])
            show(IOBuffer(), method)
            fitness = minimum( res ) 
            test_result(fitness, desired_accuracy)
            @test f_calls == res.f_calls
            @test fvals(res) == fvals(res.population)
            @test sum(abs.(hval(res.best_sol))) < 1e-3
            @test !any(gval(res.best_sol) .> 0)
            @test Metaheuristics.diff_check(res, method.information, method.options) isa Bool
            @test Metaheuristics.diversity_stop_check(res, method.information, method.options) isa Bool
        end
    end

    for problem in [:constrained1]
        run_methods_cop(problem)
    end
end
