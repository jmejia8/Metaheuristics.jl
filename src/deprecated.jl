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

function DE(f::Function, D::Int;
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


    @warn "DE(f, D;...) function is deprecated. Please use: `optimize(f, bounds, DE())`"

    a, b = limits[1,:], limits[2,:]

    if length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')

    options = Options(f_calls_limit = max_evals, store_convergence = (saveConvergence != ""))
    method = DE(;N = N,
                F = Float64(F),
               CR = Float64(CR),
           CR_min = CR_min,
            CR_max= CR_max,
            F_min =F_min,
            F_max =F_max,
        strategy  = strategy)

    method.options.f_calls_limit = max_evals

    status = optimize(f, bounds, method)

    if true
        display(status)
    end


    return status.best_sol.x, status.best_sol.f
end


function pso(f::Function, D::Int;
                        N::Int = 10D,
                       C1::Real= 2.0,
                       C2::Real= 2.0,
                        ω::Real= 0.8,
                max_evals::Int = 10000D,
              termination::Function = (x ->false),
              showResults::Bool = true,
                   limits  = (-100., 100.))

    @warn "pso(f, D;...) function is deprecated. Please use: `optimize(f, bounds, PSO())`"

    a, b = limits

    if typeof(a) <: Real
        a = ones(D) * a
        b = ones(D) * b
    elseif length(a) < D
        a = ones(D) * a[1]
        b = ones(D) * b[1]
    end

    bounds = Array([a b]')


    options = Options(f_calls_limit = max_evals)
    method = PSO(;N = N, C1 = C1, C2=C2, ω = ω)

    status = optimize(f, bounds, method)


    if showResults
        display(status)
    end

    return status.best_sol.x, status.best_sol.f
end


function ABC(
        fobj::Function,
        bounds;
        N = 50,
        limit=10,
        iters = Inf,
        max_evals = 10000*size(bounds, 2),
        termination =  x -> false,
        Ne = div(N+1, 2),
        No = div(N+1, 2),
    )
    @warn "ABC(f, bounds) is deprecated. Use optimize(f, bounds, ABC())."
    method = ABC(N = N, Ne = Ne, No = No, limit = limit)
    method.options.f_calls_limit = max_evals
    method.options.iterations = min( round(Int, max_evals / N), iters )

    method.engine.stop_criteria = (status, information, options) ->
        stop_check_abc(status, information, options) || termination(status.population)

    res = optimize(fobj, bounds, method)

    return minimizer(res), minimum(res)
end

# Gravitational Search Algorithm.
function CGSA(fobj::Function,
                D::Int;
                N::Int    = 30,
   chValueInitial::Real   = 20,
       chaosIndex::Real   = 9,
     ElitistCheck::Int    = 1,
       searchType::Symbol = :minimize,
        max_evals::Int    = 10000D,
           Rpower::Int    = 1,
         saveLast::String = "",
      showResults::Bool   = true,
      saveConvergence::String="",
              limits = [-100.0, 100.0])
	#V:   Velocity.
	#a:   Acceleration.
	#M:   Mass.  Ma = Mp = Mi = M
	#D: Dension of the test function.
	#N:   Number of agents.
	#X:   Position of agents. D-by-N matrix.
	#R:   Distance between agents in search space.
	#[low-up]: Allowable range for search space.
	#Rnorm:  Norm in eq.8.
	#Rpower: Power of R in eq.7.
	Rnorm = 2

	@warn "Deprecated function. Use optimize(f, bounds, CGSA())"

	# bounds vectors
    low, up = limits[1,:], limits[2,:]
    if length(low) < D
        low = ones(D) * low[1]
        up = ones(D) * up[1]
    end

	max_it = div(max_evals, N) + 1
	
	# random initialization for agents.
	P = initializePop(fobj, N, D, low, up)
	fitness = getfValues(P)
	X = getPositions(P, N, D)

	# Velocity
	V = zeros(N,D)

	# Current best
	theBest = getBest(P, searchType)


    convergence = []
	if saveConvergence != "" && isfeasible(theBest)
		push!(convergence, [N theBest.f])
	end

	# chaos
	wMax = chValueInitial
	wMin = 1e-10
	for iteration = 1:max_it

		# iteration
		chValue = wMax-iteration*((wMax-wMin)/max_it)
	  
		#Calculation of M. eq.14-20
		M = massCalculation(fitness,searchType)

		#Calculation of Gravitational constant. eq.13.
		G = Gconstant(iteration, max_it)

		if 1 <= chaosIndex <= 10
			G += chaos(chaosIndex,iteration,max_it,chValue)
		end

		#Calculation of accelaration in gravitational field. eq.7-10,21.
		a = Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)

		#Agent movement. eq.11-12
		X, V = move(X,a,V)

		# Checking allowable range. 
		X = correctPop(X, low, up)
		for i = 1:N
			x = X[i,:]
			P[i] = generateChild(x, fobj(x))
			fitness[i] = P[i].f
		end
		
		#Evaluation of agents. 
		currentBest = getBest(P, searchType)

		# fix this
		if Selection(theBest, currentBest, searchType)
			theBest = currentBest
		end


		if saveConvergence != "" && isfeasible(theBest)
			push!(convergence, [(iteration+1)*N theBest.f])
		end

	end #iteration

	if saveLast != ""
		writecsv(saveLast, X)        
	end

	if saveConvergence != ""
		writecsv(saveConvergence, convergence)
	end

	if showResults
		println("===========[CGSA results ]=============")
		printResults(theBest, P, max_it, max_it*N)
		println("=======================================")
	end

	return theBest.x, theBest.f

end


