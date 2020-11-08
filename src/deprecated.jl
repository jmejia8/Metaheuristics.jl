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


function SA( fobj::Function,
                D::Int;
               x0::Vector = [],
                N::Int = 500,
        max_evals::Int = 10000D,
           TolFun::Real= 1e-4,
      showResults::Bool = true,
         saveLast::String = "",
       saveConvergence::String="",
               limits  = [-100., 100.])

	@warn "Deprecated function. Use optimize(f, bounds, SA())"
	l, u = limits

	if length(x0) != D
	    x0 = l .+ (u .- l) .* rand(D)
	end

	# the current point and fx=f(x)
	x = x0
	fx= fobj(x)
	f0= fx


	nevals = 1
	stop = false

	convergence = []
	if saveConvergence != ""
		push!(convergence, [nevals f0])
	end

	t = 1
	# Main loop simulates de annealing from a high temperature to zero in max_evals.
	while !stop

		# T is the inverse of temperature.
		T = nevals / max_evals 
		μ = 10.0 ^( 100T )    
		
		# For each temperature we take 500 test points to simulate reach termal
		# equilibrium.
		for i = 1:N        
			# We generate new test point using newSol function      
			dx = newSol(2rand(length(x)) .- 1.0 , μ) .* (u-l)

			# the test point and fx1=f(x1)
			x1 = x + dx
			
			# Next step is to keep solution within bounds
			x1 = (x1 .< l).*l+(l .<= x1).*(x1 .<= u).*x1+(u .< x1).*u			
			fx1 = fobj(x1)

			nevals += 1

			df  = fx1 - fx
			
			# If the function variation,df, is <0 we take test point as current
			# point. And if df>0 we use Metropolis condition to accept or
			# reject the test point as current point.
			if (df < 0 || rand() < exp(-T*df/(abs(fx)) / TolFun))
				x = x1
				fx= fx1
			end        
			
			# If the current point is better than current solution, we take
			# current point as cuyrrent solution.       
			if fx1 < f0
				x0 = x1
				f0 = fx1
			end

			stop = nevals >= max_evals

			if stop
			    break
			end
		end

		if saveConvergence != ""
			push!(convergence, [nevals f0])
		end

		t += 1
	end

	if saveLast != ""
		writecsv(saveLast, x0)        
	end

	if saveConvergence != ""
		writecsv(saveConvergence, convergence)
	end

	if showResults
		println("===========[  SA results ]=============")
		println("| Generations = $t")
		println("| Evals       = ", nevals)
		@printf("| best f.     = %e\n", f0)
		println("=======================================")
	end

	return x0, f0
end

function WOA(fobj::Function,
                D::Int;
                N::Int = 30,
        max_evals::Int = 10000D, 
      showResults::Bool = true,
         saveLast::String = "",
  saveConvergence::String = "",
             limits = [-100.0, 100.0])

    @warn "Deprecated function. Use optimize(f, bounds, WOA())"

    Max_iter = div(max_evals, N) + 1
    lb, ub = limits

    # bounds vectors
    lb, ub = limits[1,:], limits[2,:]
    if length(lb) < D
        lb = ones(D) * lb[1]
        ub = ones(D) * ub[1]
    end


    # initialize position vector and score for the leader
    Leader_pos  = zeros(1,D)
    Leader_score= Inf #change this to -inf for maximization problems


    #Initialize the positions of search agents
    Positions = initializePop(N, D, lb, ub)

    convergence = []

    t=0# Loop counter
    nevals = 0

    # Main loop
    while t < Max_iter
        # Return back the search agents that go beyond the boundaries of the search space
        Positions = correctPop(Positions, lb, ub)

        for i=1:N
            

            # Calculate objective function for each search agent
            fitness=fobj(Positions[i,:])
            nevals += 1
            
            # Update the leader
            if fitness < Leader_score # Change this to > for maximization problem
                Leader_score = fitness # Update alpha
                Leader_pos   = Positions[i,:]
            end
            
        end
        
        a=2-t*((2)/Max_iter) # a decreases linearly fron 2 to 0 in Eq. (2.3)
        
        # a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a2=-1+t*((-1)/Max_iter)
        
        # Update the Position of search agents 
        for i=1:N
            r1 = rand() # r1 is a random number in [0,1]
            r2 = rand() # r2 is a random number in [0,1]
            
            A = 2a*r1 - a  # Eq. (2.3) in the paper
            C = 2r2        # Eq. (2.4) in the paper
            
            
            b = 1               #  parameters in Eq. (2.5)
            l = (a2-1)*rand()+1   #  parameters in Eq. (2.5)
            
            p = rand()        # p in Eq. (2.6)
            
            for j = 1:D
                
                if p < 0.5   
                    if abs(A) >= 1
                        rand_leader_index = floor(Integer, N* rand()+1)
                        X_rand  = Positions[rand_leader_index, :]
                        D_X_rand= abs(C*X_rand[j]-Positions[i,j]) # Eq. (2.7)
                        Positions[i,j]=X_rand[j]-A*D_X_rand       # Eq. (2.8)
                        
                    elseif abs(A)<1
                        D_Leader=abs(C*Leader_pos[j]-Positions[i,j]) # Eq. (2.1)
                        Positions[i,j]=Leader_pos[j]-A*D_Leader      # Eq. (2.2)
                    end
                    
                elseif p>=0.5
                  
                    distance2Leader=abs(Leader_pos[j]-Positions[i,j])
                    # Eq. (2.5)
                    Positions[i,j]=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos[j]
                    
                end
                
            end
        end
        t=t+1

        push!(convergence, [nevals Leader_score])

    end

    if saveConvergence != ""
       writecsv(saveConvergence, convergence)
    end

    if saveLast != ""
       writecsv(saveLast, Positions)
    end
    
    if showResults
        println("===========[ ECA results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", t*N)
        @printf("| best f.     = %e\n", Leader_score)
        println("=======================================")
    end

    return Leader_pos, Leader_score


end


