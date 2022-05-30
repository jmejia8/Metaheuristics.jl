using JMcDM

@testset "DecisionMaking: JMcDM" begin
    function test_jmcdm()
        _, _, pf = Metaheuristics.TestProblems.ZDT1();
        res = State(pf[1], pf)

        w = fill(0.5, 2);

        # PrometheeMethod
        # PSIMethod, MoosraMethod --err
        methods = [
                   ArasMethod ,CocosoMethod ,CodasMethod ,CoprasMethod ,CriticMethod,
                   EdasMethod ,ElectreMethod ,GreyMethod ,MabacMethod ,MaircaMethod ,MooraMethod,
                   SawMethod ,TopsisMethod ,VikorMethod ,WPMMethod ,WaspasMethod,
                   MarcosMethod
                  ]

        for method in methods
            res_dm = mcdm(res, w, method())
            res_dm2 = mcdm(MCDMSetting(res, w), method())
            @test res_dm.bestIndex == res_dm2.bestIndex
        end


        # JMcDM.summary(res, w, methods)
    end

    test_jmcdm()

end


