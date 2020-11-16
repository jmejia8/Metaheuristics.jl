# Gravitational Search Algorithm.
function GSA(fobj::Function,
                D::Int;
                N::Int    = 30,
     ElitistCheck::Int    = 1,
       searchType::Symbol = :minimize,
        max_evals::Int    = 10000D,
           Rpower::Int    = 1,
         saveLast::String = "",
      showResults::Bool = true,
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

	low, up = limits

	# bounds vectors
	low = low * ones(D)
	up  =  up * ones(D)

	max_it = div(max_evals, N) + 1

	#random initialization for agents.
	X = initializePop(N, D, low, up)

	#Evaluation of agents. 
	fitness = evaluatePop(X, fobj, N)

	V = zeros(N,D)

	# Current best
	best_X, best = getBest(fitness, searchType)

	Fbest = best
	Lbest = X[best_X,:]

    convergence = []
	if saveConvergence != ""
		push!(convergence, [N Fbest])
	end

	for iteration = 1:max_it
		# iteration
	  
		#Calculation of M. eq.14-20
		M = massCalculation(fitness,searchType)

		#Calculation of Gravitational constant. eq.13.
		G = Gconstant(iteration, max_it)

		#Calculation of accelaration in gravitational field. eq.7-10,21.
		a = Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)

		#Agent movement. eq.11-12
		X, V = move(X,a,V)

		# Checking allowable range. 
		X = correctPop(X, low, up)
		
		#Evaluation of agents. 
		fitness = evaluatePop(X, fobj, N)
		best_X, best = getBest(fitness, searchType)


		if searchType == :minimize
			if best < Fbest  # minimization.
				Fbest = best
				Lbest = X[best_X,:]
			end
		else 
			if best > Fbest  # maximization
				Fbest = best
				Lbest = X[best_X,:]
			end
		end

		if saveConvergence != ""
			push!(convergence, [(iteration+1)*N Fbest])
		end

	end #iteration

	if saveLast != ""
		writecsv(saveLast, X)        
	end

	if saveConvergence != ""
		writecsv(saveConvergence, convergence)
	end

	if showResults
		println("===========[ GSA results ]=============")
		println("| Generations = $max_it")
		println("| Evals       = ", max_it*N)
		@printf("| best f.     = %e\n", Fbest)
		println("=======================================")
	end

	return Lbest, Fbest

end
