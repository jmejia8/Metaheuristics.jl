function eca(fobj::Function,
                D::Int;
            η_max::Real  = 2,
                K::Int   = 7,
                N::Int   = K*D,
        p_exploit::Real  = 0.95,
            p_bin::Real  = 0.02,
        max_evals::Int   = 10000D,
      showResults::Bool  = true,
       correctSol::Bool  = true,
       searchType::Symbol=:minimize,
      initPopRand::Symbol=:uniform,
         showIter::Bool  = false,
         saveLast::String= "",
         adaptive::Bool  = false,
     canResizePop::Bool  = false,
      termination::Function   = (x ->false),
       saveConvergence::String="",
    returnDetails::Bool = false,
           limits  = [-100., 100.])

    @warn "eca(f, D;...) function is deprecated. Please use: `optimize(f, bounds, ECA())`"

    f = fobj
    
    a, b = limits[1,:], limits[2,:]

    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')

    options = Options(f_calls_limit = max_evals, debug=showIter, store_convergence = (saveConvergence != ""))
    method = ECA()

    method.options.f_calls_limit = max_evals

    status = optimize(f, bounds, method)

    if showResults
        display(status)
    end

    return status.best_sol.x, status.best_sol.f
end

function diffEvolution(func::Function, D::Int;
                        N::Int = 10D,
                        F::Real= 1.0,
                       CR::Real= 0.9,
                   CR_min::Real= CR,
                   CR_max::Real= CR,
                    F_min::Real=F,
                    F_max::Real=F,
                max_evals::Int = 10000D,
                 strategy::Symbol = :rand1,
              termination::Function = (x ->false),
              showResults::Bool = true,
                 saveLast::String = "",
          saveConvergence::String="",
                   limits  = [-100., 100.])
    warn("This function is deprecated. Use DE function")
    return DE(func, D;
                    N  = N,
                    F  = F,
                    CR = CR,
                    CR_min = CR_min,
                    CR_max = CR_max,
                    F_min = F_min,
                    F_max = F_max,
                    max_evals  = max_evals,
                    strategy  = strategy,
                    termination  = termination,
                    showResults  = showResults,
                    saveLast  = saveLast,
                    saveConvergence = saveConvergence,
                    limits  = limits)

end
