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

        # handling empty constraints
        sol = Metaheuristics.create_child(x, (1.0, zeros(0), zeros(2)))
        @test !isnan(Metaheuristics.sum_violations(sol)) && Metaheuristics.sum_violations(sol) == 0
        sol = Metaheuristics.create_child(x, (1.0, -ones(1), zeros(0)))
        @test !isnan(Metaheuristics.sum_violations(sol)) && Metaheuristics.sum_violations(sol) == 0
        sol = Metaheuristics.create_child(x, (1.0, zeros(0), zeros(0)))
        @test sol isa Metaheuristics.xf_indiv
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
        hyperVolume = Metaheuristics.PerformanceIndicators.hypervolume(front, referencePoint, verbose=false)
        @test hyperVolume ≈ 6.793879034744429

        front_ = Array(hcat(front...)')
        hyperVolume = Metaheuristics.PerformanceIndicators.hypervolume(front_, referencePoint)
        @test hyperVolume ≈ 6.793879034744429

        # testing nadir and ideal pint
        @test ideal(front) == ideal(Array(hcat(front...)')) == zeros(3)
        @test nadir(front) == nadir(Array(hcat(front...)')) == [2.0, 3, 2]
        @test isempty(nadir(empty(front))) && isempty(ideal(empty(front)))

        f, bounds, front2 = Metaheuristics.TestProblems.get_problem(:DTLZ2);
        @test Metaheuristics.PerformanceIndicators.hypervolume(front2, nadir(front2)) > 0
        @test ideal(front2) != nadir(front2) 
        @test isempty(nadir(empty(front2))) && isempty(ideal(empty(front2)))

        for  s in front2
            s.sum_violations = 1
        end

        # ideal and nadir points on an infeasible front
        @test  Metaheuristics.compare(ideal(front2), nadir(front2)) == 1
        
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

    function test_reproduction()
        D = 4
        N = 10
        # Parallel Evaluations of solutions
        pop1 = Metaheuristics.create_child(rand(N,D), rand(N)) 
        pop2 = Metaheuristics.create_child(rand(N,D), (rand(N), -ones(N,2), zeros(N,2))) 

        problem = Metaheuristics.Problem(x -> rand(), [zeros(D)'; ones(D)'])
        for pop in [pop1, pop2]
            status = State(pop[1], pop)
            X = Metaheuristics.reproduction(status, ECA(N=N,K=3).parameters, problem)
            @test size(X) == (N, D)
            @test length(pop) == N
        end

        v = Metaheuristics.GA_reproduction_half(rand(D), rand(D), problem.bounds)
        @test length(v) == D
    end

    function test_sampling_methods()
        seed!(1)
        # current methods implemented for sampling
        methods = [LatinHypercubeSampling(100, 2), Grid(10, 2), RandomInBounds(100)]

        # bounds
        a = [-10 -10.]
        b = [10 10.0]
        bounds = [a; b]
        # perform tests
        for (i, method) in enumerate(methods)
            if !(method isa RandomInBounds) # RandomInBounds required bounds
                # test if values are in [0,1]
                @test !any( .!(0 .<= sample(method) .<= 1))
            end
            X = sample(method, bounds)
            # check the minimum pairwise distance
            @test Metaheuristics._score_lhs(X) > 0
            # checking if points are out of bounds
            @test !any( .!(a .<= X .<= b) )
        end
    end

    function test_performance_indicators()
        A1 = [4 7;5 6;7 5; 8 4.0; 9 2]
        A2 = [4 7;5 6;7 5; 8 4.0]
        A3 = [6 8; 7 7;8 6; 9 5;10 4. ]

        @test PerformanceIndicators.epsilon_indicator(A1, A2) == 1
        @test PerformanceIndicators.epsilon_indicator(A1, A3) == 9/10
        @test PerformanceIndicators.epsilon_indicator(A3, A2) == 3/2

        # check for translated fronts
        @test PerformanceIndicators.epsilon_indicator(A1 .- 10, A1 .- 10) == 1

        _, _, pf = Metaheuristics.TestProblems.ZDT3()
        @test PerformanceIndicators.epsilon_indicator(pf, pf) == 1
    end
    

    test_reproduction()
    simple_test()
    test_problems()
    test_prev_pop()
    test_hypervolume()
    test_sampling_methods()
    test_performance_indicators()

end
