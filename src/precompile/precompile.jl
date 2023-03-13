using SnoopPrecompile

@precompile_setup begin



    @precompile_all_calls begin
        
        begin
            local f(x) = sum(x .^ 2)
            optimize(f, [-1 -1 -1; 1 1 1.0], ABC())
            optimize(f, [-1 -1 -1; 1 1 1.0], CGSA())
            optimize(f, [-1 -1 -1; 1 1 1.0], DE())
            optimize(f, [-1 -1 -1; 1 1 1.0], ECA())
            optimize(f, [-1 -1 -1; 1 1 1.0], PSO())
            optimize(f, [-1 -1 -1; 1 1 1.0], SA())
            optimize(f, [-1 -1 -1; 1 1 1.0], WOA())
        end

        begin
            local f, bounds, pf = Metaheuristics.TestProblems.MTP()
            ccmo = CCMO(NSGA2(N = 100, p_m = 0.001))
            optimize(f, bounds, ccmo)
        end

        begin
            local g(x) = sum(abs.(x .- (length(x):-1:1.0)))
            perm_size = 10
            ga = GA(;
                initializer = RandomPermutation(N = 100),
                crossover = OrderCrossover(),
                mutation = SlightMutation(),
            )
            optimize(g, zeros(Int, 2, perm_size), ga)
        end

        begin
            local f, bounds, solutions = Metaheuristics.TestProblems.rastrigin()
            result = optimize(f, bounds, MCCGA())
        end

        begin
            D = 2
            local h(x) = (x, [sum(x .^ 2) - 1], [0.0])
            bounds = [
                -1 -1
                1 1.0
            ]
            nobjectives = 2
            npartitions = 100
            weights = gen_ref_dirs(nobjectives, npartitions)
            moead_de = MOEAD_DE(weights, options = Options(debug = false, iterations = 250))
            status_moead = optimize(h, bounds, moead_de)
        end

        begin
            D = 2
            local j(x) = (x, [sum(x .^ 2) - 1], [0.0])
            bounds = [
                -1 -1
                1 1.0
            ]
            nsga2 = NSGA2(N = 100, p_cr = 0.85)
            status = optimize(j, bounds, nsga2)
        end


        begin
            local g, bounds, pf = Metaheuristics.TestProblems.get_problem(:DTLZ2)
            nsga3 = NSGA3(p_cr = 0.9)
            status = optimize(g, bounds, nsga3)
        end

        begin
            local k(x) = (x, [sum(x .^ 2) - 1], [0.0])
            bounds = [
                -1 -1
                1 1.0
            ]
            sms_emoa = SMS_EMOA(N = 100, p_cr = 0.85)
            status = optimize(k, bounds, sms_emoa)
        end

        begin
            local m(x) = (x, [sum(x .^ 2) - 1], [0.0])
            bounds = [
                -1 -1
                1 1.0
            ]
            nsga2 = SPEA2(N = 100, p_cr = 0.85)
            status = optimize(m, bounds, nsga2)
        end

    end
end
