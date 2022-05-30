using JMcDM

@testset "DecisionMaking: ROI" begin
    function test_roi()
        _, _, pf = Metaheuristics.TestProblems.ZDT1();
        res = State(pf[1], pf)

        w = [0.1 0.9; 0.9 0.1]
        δ = [0.1,0.1]
        method = ROIArchiving(δ)
        idx = decisionmaking(res, w, method)
        @test !isempty(idx)
        sols = best_alternative(res, w, method)
        @test length(sols) == length(idx)
    end
    
    test_roi() 
end

@testset "DecisionMaking: JMcDM" begin
    function test_jmcdm()
        _, _, pf = Metaheuristics.TestProblems.ZDT1();
        res = State(pf[1], pf)

        w = [0.5, 0.5]

        # PrometheeMethod
        # PSIMethod, MoosraMethod --err
        methods = [
                   ArasMethod, CocosoMethod, CodasMethod, CoprasMethod, 
                   EdasMethod, ElectreMethod, GreyMethod, MabacMethod, MaircaMethod,
                   MooraMethod, SawMethod, TopsisMethod, VikorMethod, WPMMethod,
                   WaspasMethod, MarcosMethod
                  ]

        for method in methods
            res_dm = mcdm(res, w, method())
            res_dm2 = mcdm(MCDMSetting(res, w), method())
            @test res_dm.bestIndex == res_dm2.bestIndex
            
            # bestIndex can be touple and needs to be handled...
            idx = res_dm.bestIndex isa Tuple ? res_dm.bestIndex[1] : res_dm.bestIndex
            ref_sol = res.population[idx]
            best_sol_ = best_alternative(res, w, method())
            best_sol = best_sol_ isa Array ? best_sol_[1] : best_sol_
            @test Metaheuristics.compare(fval(best_sol), fval(ref_sol)) == 0
        end

        @test JMcDM.summary(res, w, [:topsis, :electre, :vikor]) isa DataFrame

    end

    test_jmcdm()

end

