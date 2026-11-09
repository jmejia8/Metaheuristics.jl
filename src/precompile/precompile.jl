using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    opts = Options(iterations = 4, f_calls_limit = 1000, seed = 1)
    fso(x) = sum(abs2, x)
    fco(x) = (sum(abs2, x), [x[1] - 0.5], [0.0])
    fmo(x) = ([x[1], 1 - sqrt(abs(x[1])) + sum(x[2:end])], [0.0], [0.0])
    fh(x) = (x, [sum(x.^2) - 1], [0.0])
    bmat = [-1 -1 -1; 1 1 1.0]
    bbox = boxconstraints(lb = -ones(3), ub = ones(3))
    bm01 = boxconstraints(lb = zeros(3), ub = ones(3))
    @compile_workload begin
        for b in (bmat, bbox), mk in (ABC, CGSA, DE, ECA, PSO, SA, WOA, SHADE, CSO)
            opts = Options(iterations = 4, f_calls_limit = 1000, seed = 1)
            optimize(fso, b, mk(N = 10, options = opts))
        end
        for b in (bmat, bbox)
            optimize(fso, b, RDEx(options = opts))
        end
        for b in (bmat, bbox), mk in (ECA, DE, εDE)
            r = optimize(fco, b, mk(N = 10, options = opts))
            sprint(show, MIME"text/plain"(), r)
        end
        r = optimize(fso, bbox, ECA(N = 10, options = opts))
        sprint(show, MIME"text/plain"(), r)
        optimize(fso, bm01, BRKGA(options = opts))
        let (f, bounds, pf) = Metaheuristics.TestProblems.MTP()
            optimize(f, bounds, CCMO(NSGA2(N = 100, p_m = 0.001), options = opts))
        end
        let g(x) = sum(abs.(x .- (length(x):-1:1.0)))
            optimize(
                g,
                zeros(Int, 2, 10),
                GA(;
                    initializer = RandomPermutation(N = 100),
                    crossover = OrderCrossover(),
                    mutation = SlightMutation(),
                    options = opts,
                ),
            )
        end
        let (f, bounds, _) = Metaheuristics.TestProblems.rastrigin()
            r = optimize(f, bounds, MCCGA(options = opts))
            fvals(r)
            positions(r)
            minimizer(r)
            minimum(r)
        end
        for b in ([-1 -1; 1 1.0], boxconstraints(lb = -ones(2), ub = ones(2)))
            optimize(
                fh,
                b,
                MOEAD_DE(gen_ref_dirs(2, 100), options = Options(iterations = 2)),
            )
        end
        for b in ([-1 -1; 1 1.0], bm01), mk in (NSGA2, SMS_EMOA, SPEA2)
            r = optimize(
                b isa Matrix ? fh : fmo,
                b,
                mk(N = 100, p_cr = 0.85, options = opts),
            )
            fvals(r)
            sprint(show, MIME"text/plain"(), r)
            PerformanceIndicators.hypervolume(pareto_front(r), [2.0, 2.0])
        end
        let (g, bounds, _) = Metaheuristics.TestProblems.get_problem(:DTLZ2)
            optimize(g, bounds, NSGA3(p_cr = 0.9, options = opts))
        end
    end
end
