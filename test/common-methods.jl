using Metaheuristics
using Test

@testset "Common Methods" begin
    function simple_test()

        # compare
        x = rand(10)

        @test Metaheuristics.compare(x, x) == 0
        @test Metaheuristics.compare(x, x .+ 1) == 1
        @test Metaheuristics.compare(x .+ 1, x) == 2
        y = copy(x)
        y[1] -= 1.0 
        x[end] -= 1.0
        @test Metaheuristics.compare(x, y) == 3

    end 

    function test_problems()

        # single-objective
        for problem = [:sphere, :discus, :rastrigin, ]
            f, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)
            for opt in optimums
                @test f(opt.x) == opt.f == 0.0
            end
        end


        # single-objective constrained
        for problem = [:constrained1, :constrained2, :constrained3]
            f, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)
            for opt in optimums
                fx, gx, hx = f(opt.x)
                @test fx == opt.f && sum(max.(gx, 0.0)) ≈ 0.0 && sum(abs.(hx)) ≈ 0.0
            end
        end
        
        # multi-objective
        for problem = [:ZDT1, :ZDT2, :ZDT3, :ZDT4, :ZDT6, :DTLZ2]
            f, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)
            pf = Metaheuristics.get_non_dominated_solutions(optimums)

            @test length(pf) == length(optimums)
            
            # test dimensions
            fx, gx, hx = f(bounds[1,:])
            @test length(fx) == length(optimums[1].f)

            # test performance indicators
            @test Metaheuristics.PerformanceIndicators.igd(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.igd_plus(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.gd(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.gd_plus(optimums, pf) ≈ 0.0
        end

	
    end
    

    simple_test()
    test_problems()
end
