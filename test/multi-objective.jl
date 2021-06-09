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
            f_calls += 1
            ff(x)
        end

        options = Options( seed = 1, iterations=10000, f_calls_limit = 25000)
        nobjectives = length(pf[1].f)
        npartitions = nobjectives == 2 ? 100 : 12

        methods = [
                NSGA2(options=options),
                MOEAD_DE(gen_ref_dirs(nobjectives, npartitions), options=Options( seed = 1, iterations = 500)),
                NSGA3(options=options),
              ]

        for method in methods
            f_calls = 0
            result = ( optimize(f, bounds, method) )
            show(IOBuffer(), "text/html", result)
            show(IOBuffer(), "text/plain", result.population)
            show(IOBuffer(), result.population[1])
            @test Metaheuristics.PerformanceIndicators.igd(result.population, pf) <= 0.2
            @test Metaheuristics.PerformanceIndicators.spacing(result) < 0.2
            @test Metaheuristics.PerformanceIndicators.covering(pf, result.population) <= 1.0
            @test Metaheuristics.PerformanceIndicators.covering(result, result) ≈ 0.0

            # number of function evaluations should be reported correctly
            @test f_calls == result.f_calls


            # test obtaining non-dominated solutions
            pf1 = pareto_front(result)
            pf2 = pareto_front(result.population)

            @test size(pf1, 1) == size(pf2,1) &&
                  Metaheuristics.PerformanceIndicators.igd(pf1, pf2) ≈ 0.0
        end
    end


    for problem in [:ZDT3, :DTLZ2, :MTP]
        run_methods(problem)
    end




end
