@testset "Initial Solutions" begin

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
        x_optimum = get_position(first(optimums))

        methods = [ABC, DE, ECA, PSO, WOA, MCCGA, CGSA]

        for method in methods
            for N in [5, 10, 20] # vary the population size
                algo = method(N = N, options = options, information = information)
                # initial solutions
                X = rand(10, size(bounds, 2))
                X[1,:] = x_optimum #zeros(size(bounds, 2))
                set_user_solutions!(algo, X, f)

                res = optimize(f, bounds, algo)

                @test length(res.population) == algo.parameters.N 
                @test minimum(res) ≈ first(fvals(optimums))
                @test sum(abs,minimizer(res)-x_optimum) ≈ 0
                # TODO check the number of function evaluations
                #@test f_calls == nfes(res)
            end
        end
    end

    function run_initialized_multiobjective(problem)
        ff, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        options = Options(iterations = 3, f_calls_limit = 500)
        information = Information()
        x_optimum = get_position(first(optimums))

        methods = [SPEA2, SMS_EMOA, NSGA2, NSGA3]

        for method in methods
            for N in [5, 10, 20] # vary the population size
                algo = method(N = N, options = options, information = information)
                # initial solutions
                X = rand(10, size(bounds, 2))
                X[1,:] = x_optimum
                set_user_solutions!(algo, X, f)

                res = optimize(f, bounds, algo)
                @test length(res.population) == algo.parameters.N 
                # check whether optimum was used by optimizer
                v = sum(abs.(positions(res) .- x_optimum'), dims=2) |> minimum 
                @test  v ≈ 0
            end
        end
    end

    for problem in [:sphere]
        run_initialized(problem)
    end

    for problem in [:ZDT1]
        run_initialized_multiobjective(problem)
    end
end
