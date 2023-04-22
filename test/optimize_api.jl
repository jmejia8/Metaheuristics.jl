@testset "Optimize API" begin
    function test_optimize()
        f, bounds, _ = Metaheuristics.TestProblems.sphere();
        method = ECA()
        while !Metaheuristics.should_stop(method)
            optimize!(f, bounds, method)
        end
        result = Metaheuristics.get_result(method)
        @test !isnothing(result.best_sol)
    end

    function search_space_optimize()
        f, _, _ = Metaheuristics.TestProblems.sphere();
        lb = zeros(3)
        ub = ones(3)
        options = Options(iterations = 10, seed = 1)
        show(IOBuffer(), "text/plain", options)
        def_bounds = [
                         (lb, ub),
                         boxconstraints(lb, ub, rigid=false),
                         [lb ub],
                         [lb ub]',
                         [-5; 5;;], # due to #83
                        ]
                        
        for space in def_bounds
            res = optimize(f, space, ECA(;K = 3, options))
            res2 = optimize(f, space, ECA, K = 3, iterations = 10, seed = 1)
            @test minimum(res) == minimum(res2)
        end
        for space in def_bounds
            p = Metaheuristics.Problem(sum, space)
            for algo in [ABC, DE, PSO, WOA, CGSA, GA, SA]
                res = optimize(f, space, algo, iterations=3, time_limit=0.1)
                x = minimizer(res)
                @test x isa AbstractVector
                @test length(x) == Metaheuristics.getdim(p)
            end
            for algo in [SPEA2, SMS_EMOA, NSGA2, NSGA3]
                res = optimize(x->([0.0,0],[0.0],[0.0]),
                               space, algo, iterations=3, time_limit=0.1
                              )
                @test size(positions(res), 2) == Metaheuristics.getdim(p)
                @test eltype(positions(res)) <: eltype(p.search_space.lb)
            end
        end

        redirect_stdout(devnull) do
            optimize(f, [lb ub], ECA, iterations=10, verbose = true)
        end
    end
    test_optimize()
    search_space_optimize()
end
