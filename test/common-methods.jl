using Metaheuristics
using Test
import Random: seed!

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

        sol_feasible   = Metaheuristics.create_child(x, (x, zeros(2), zeros(2)))
        sol_infeasible = Metaheuristics.create_child(x, (x, rand(4), rand(8)))

        @test Metaheuristics.compare(sol_feasible, sol_infeasible) == 1
        @test Metaheuristics.compare(sol_infeasible, sol_feasible) == 2
        @test Metaheuristics.is_better(sol_feasible, sol_infeasible) # true
        @test !Metaheuristics.is_feasible(sol_infeasible) # false
        @test Metaheuristics.is_feasible(sol_feasible) # true
        @test !Metaheuristics.is_better(sol_infeasible, sol_feasible)# !false == true

        # repair
        x = collect(range(-10, 10, length=20))
        bounds = [-ones(20) ones(20)]'
        for repair in [Metaheuristics.reflected_back_to_bound!,
                       Metaheuristics.replace_with_random_in_bounds!,
                       Metaheuristics.wrap_to_bounds!,
                       Metaheuristics.reset_to_violated_bounds!,
                      ]
            @test Metaheuristics.is_in_bounds(repair(x, bounds), bounds)
        end
        @test Metaheuristics.is_in_bounds(Metaheuristics.evo_boundary_repairer!(x, zeros(20), bounds), bounds)

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
        for problem = [:ZDT1, :ZDT2, :ZDT3, :ZDT4, :ZDT6,
                       :DTLZ1, :DTLZ2, :DTLZ3, :DTLZ4, :DTLZ5, :DTLZ6,
                       :C1_DTLZ1, :C1_DTLZ3, :C2_DTLZ2, :C3_DTLZ4
                      ]
            f, bounds, optimums = Metaheuristics.TestProblems.get_problem(problem)
            pf = Metaheuristics.get_non_dominated_solutions(optimums)

            @test length(pf) == length(optimums)
            
            # test dimensions
            fx, gx, hx = f(bounds[1,:])
            @test length(fx) == length(optimums[1].f)

            # test performance indicators
            ## generational distance
            @test Metaheuristics.PerformanceIndicators.igd(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.igd_plus(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.gd(optimums, pf) ≈ 0.0
            @test Metaheuristics.PerformanceIndicators.gd_plus(optimums, pf) ≈ 0.0


            ## spacing
            seed!(1)
            X = rand(100, 7)

            @test Metaheuristics.PerformanceIndicators.spacing(X) < 0.5
            @test Metaheuristics.PerformanceIndicators.spacing(optimums) < 0.2

            ## covering
            ## non dominated solutions
            X = fvals(pf)
            Y = copy(X)
            # 50% of solutions in Y will be dominated by X
            mask = 1:size(X,1) ÷ 2
            Y[mask,1] = Y[mask,1] .+ 1.0
            Cxy=Metaheuristics.PerformanceIndicators.covering(X, Y)
            Cyx=Metaheuristics.PerformanceIndicators.covering(Y, X)
            @test Cxy == (size(X,1) ÷ 2) / size(X,1)
            @test Cyx == 0.0
 
            ## Δₚ (delta p)
            @test PerformanceIndicators.deltap(pf, optimums) == PerformanceIndicators.Δₚ(pf, optimums)
        end

	
    end


    function test_hypervolume()
        front = [[0.        ,0.        ,1.1       ],
                 [0.        ,1.07863874,0.21572775],
                 [0.13540064,0.94780448,0.54160256],
                 [0.29938208,0.7484552 ,0.7484552 ],
                 [0.46669048,0.62225397,0.77781746],
                 [0.6350853 ,0.6350853 ,0.6350853 ],
                 [0.7484552 ,0.7484552 ,0.29938208],
                 [0.89510682,0.        ,0.63936201],
                 [0.98386991,0.49193496,0.        ],
                 [2.0, 3, 2], # only for testing
                 [1.1       ,0.        ,0.        ]]

        referencePoint = [2.0, 2, 2]
        hyperVolume = Metaheuristics.PerformanceIndicators.hypervolume(front, referencePoint)
        @test hyperVolume ≈ 6.793879034744429

        front_ = Array(hcat(front...)')
        hyperVolume = Metaheuristics.PerformanceIndicators.hypervolume(front_, referencePoint)
        @test hyperVolume ≈ 6.793879034744429

        # testing nadir and ideal pint
        @test ideal(front) == ideal(Array(hcat(front...)')) == zeros(3)
        @test nadir(front) == nadir(Array(hcat(front...)')) == [2.0, 3, 2]

        f, bounds, front2 = Metaheuristics.TestProblems.get_problem(:DTLZ2);
        @test Metaheuristics.PerformanceIndicators.hypervolume(front2, nadir(front2)) > 0
        @test ideal(front2) != nadir(front2) 
    end

    function test_prev_pop()
        f, bounds, optimums = Metaheuristics.TestProblems.get_problem(:sphere)

        options=Options(f_calls_limit=10000, seed=1)

        method = DE(options = options)
        res = optimize(f, bounds, method)
        optimum = minimum(res)
        N = length(res.population)

        methods = [ECA(N=N,options = options),
                   PSO(N=N,options = options),
                   DE(N=N,options = options)]

        for method in methods
            prev_status = State(res.best_sol, res.population)
            method.status = prev_status
            res = optimize(f, bounds, method)
            optimum = minimum(res)
        end
        
        @test optimum < 1e-5
    end
    
    

    simple_test()
    test_problems()
    test_prev_pop()
    test_hypervolume()

end
