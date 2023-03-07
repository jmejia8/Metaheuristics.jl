
@testset "Combinatorial" begin

    function run_initialized(problem)
        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        options = Options(iterations = 3, f_calls_limit = 500)
        information = Information()

        methods = [ABC, DE, ECA, PSO, SA, WOA, MCCGA, CGSA]

        for method in methods
            println(method)
            for N in [5, 10, 20] # vary the population size
                println(N)
                algo = method(N = N, options = options, information = information)

                # initial solutions
                X = rand(10, size(bounds, 2))
                set_inital_solutions!(algo, X, f)

                res = optimize(f, bounds, algo)
                @test length(res.population) == N #max(N, size(X, 1))
                #@test f_calls == nfes(res)
            end
            println("-----")
        end
    end

    for problem in [:sphere]
        run_initialized(problem)
    end
end
