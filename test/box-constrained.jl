using Metaheuristics
using Test



@testset "Box-constrained" begin


    function run_methods(problem)

        desired_accuracy = 1e-4

        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        information = Information(f_optimum = 0.0)
        options = Options(f_tol = desired_accuracy, seed = 1)
        options_2 = Options(f_tol = desired_accuracy, seed = 1, store_convergence=true)

        methods = [
                   ABC(options = options, information = information), 

                   # differential evolution
                   DE(CR = 0.5, options = options_2, information = information),
                   DE(CR = 0.4, strategy=:best1, options = options_2, information = information),
                   DE(CR = 0.3, strategy=:rand2, options = options_2, information = information),
                   DE(CR = 0.5, strategy=:randToBest1, options = options_2, information = information),
                   DE(CR = 0.3, strategy=:best2, options = options_2, information = information),

                   ECA(options = options, information = information),
                   ECA(adaptive=true, resize_population=true, options = options, information = information),

                   PSO(options = options, information = information),
                   SA(options = options, information = information),
                   WOA(options = options, information = information),
                  ]

        # CGSA
        for i = 1:10
            push!(methods, CGSA(chaosIndex=i,options = options, information = information))
        end

        for method in methods
            f_calls = 0
            res = optimize(f, bounds, method)
            show(IOBuffer(), res)
            show(IOBuffer(), "text/plain", res.population)
            show(IOBuffer(), "text/html", res.population)
            show(IOBuffer(), res.population[1])
            fitness = minimum( res ) 
            test_result(fitness, desired_accuracy)

            # also count the mumber of function evaluations
            @test f_calls == nfes(res)
        end
    end


    function run_methods_parallel_eval(problem)

        desired_accuracy = 1e-4

        f_calls = 0
        f(X) = begin
            f_calls += size(X,1)
            sum(X.^2,dims=2)[:,1]
        end
        bounds = [-ones(5)'; ones(5)']
        options = Options(f_tol = desired_accuracy,
                          f_calls_limit=1000,
                          seed = 1,
                          parallel_evaluation=true
                         )

        methods = [
                   CGSA(options = options),
                   DE(options = options),
                   ECA(options = options),
                   PSO(options = options),
                   WOA(options = options),
                  ]
    
        for method in methods
            f_calls = 0
            res = optimize(f, bounds, method)
            @test f_calls == nfes(res)
        end
    end

    function test_result(fitness, tol)
        @test ≈(fitness, 0.0, atol=tol)
    end

    for problem in [:sphere]
        run_methods(problem)
        run_methods_parallel_eval(problem)
    end

end
