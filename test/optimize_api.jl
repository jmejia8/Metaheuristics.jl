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
        for space in [ (lb, ub), boxconstraints(lb, ub, rigid=false), [lb ub], [lb ub]']
            res = optimize(f, space, ECA(;K = 3, options))
            res2 = optimize(f, space, ECA, K = 3, iterations = 10, seed = 1)
            @test minimum(res) == minimum(res2)
        end

        redirect_stdout(devnull) do
            optimize(f, [lb ub], ECA, iterations=10, verbose = true)
        end
    end
    test_optimize()
    search_space_optimize()
end
