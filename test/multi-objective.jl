using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "Multi objective" begin

    function run_methods(problem)


        f, bounds, pf = Metaheuristics.TestProblems.get_problem(problem)
        D = size(bounds, 2)

        options = Options( seed = 1, iterations = 500)

        methods = [
                ECA(N = 100, options=options),
                MOEAD_DE(D, 2, N = 300, options=options),
                NSGA2(options=options)
              ]

        for method in methods
            result = ( optimize(f, bounds, method) ) 
            @show Metaheuristics.PerformanceIndicators.igd(fvals(result), fvals(pf))
            display(result)
            @test true
        end
    end


    for problem in [:ZDT3]
        run_methods(problem) 
    end




end
