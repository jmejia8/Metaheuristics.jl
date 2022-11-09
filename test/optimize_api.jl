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
    test_optimize()
end
