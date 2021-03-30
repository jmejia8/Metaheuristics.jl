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


        ff, bounds, pf = Metaheuristics.TestProblems.get_problem(problem)
        D = size(bounds, 2)

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            global f_calls += 1
            ff(x)
        end

        options = Options( seed = 1, iterations = 200)

        methods = [
                # ECA(N = 100, options=options),
                MOEAD_DE(D, 2, N = 300, options=options),
                NSGA2(options=options)
              ]

        for method in methods
            global f_calls = 0
            result = ( optimize(f, bounds, method) ) 
            @test Metaheuristics.PerformanceIndicators.igd(result.population, pf) <= 0.2
            @test f_calls == result.f_calls
        end
    end


    for problem in [:ZDT3]
        run_methods(problem) 
    end




end
