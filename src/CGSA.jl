# GSA code v0.1.
# Coded by Jesús Mejía. 
# Based on MATLAB code of Esmat Rashedi, 2010. 
# Adopted from 
# https://la.mathworks.com/matlabcentral/fileexchange/61116-gsa-+-chaotic-gravitational-constant
# " E. Rashedi, H. Nezamabadi-pour and S. Saryazdi,
# “GSA: A Gravitational Search Algorithm”, Information sciences, vol. 179,
# no. 13, pp. 2232-2248, 2009."
# ------------------------------------------
# Mirjalili, Seyedali, and Amir H. Gandomi. 
# "Chaotic gravitational constants for the gravitational search algorithm." 
# Applied Soft Computing 53 (2017): 407-419.
# ------------------------------------------
# Common functions
# include("tools.jl")
# include("CGSA.jl")

function chaos(index,curr_iter,max_iter,Value)
	x = zeros(max_iter + 1)
	x[1]=0.7

	G = zeros(max_iter)
	if index == 1
		# Chebyshev map
		for i=1:max_iter
			x[i+1]=cos(i*acos(x[i]))
			G[i]=((x[i]+1)*Value)/2
		end
	elseif index == 2
		# Circle map
		a=0.5
		b=0.2
		for i=1:max_iter
			x[i+1]=mod(x[i]+b-(a/(2*pi))*sin(2*pi*x[i]),1)
			G[i]=x[i]*Value
		end
	elseif index == 3
		# Gauss/mouse map
		for i=1:max_iter
			if x[i]==0
				x[i+1]=0
			else
				x[i+1]=mod(1/x[i],1)
			end
			G[i]=x[i]*Value
		end
	elseif index == 4
		# Iterative map
		a=0.7
		for i=1:max_iter
			x[i+1]=sin((a*pi)/x[i])
			G[i]=((x[i]+1)*Value)/2
		end
	elseif index == 5
		# Logistic map
		a=4
		for i=1:max_iter
			x[i+1]=a*x[i]*(1-x[i])
			G[i]=x[i]*Value
		end
	elseif index == 6
		# Piecewise map
		P=0.4
		for i=1:max_iter
			if x[i]>=0 && x[i]<P
				x[i+1]=x[i]/P
			end
			if x[i]>=P && x[i]<0.5
				x[i+1]=(x[i]-P)/(0.5-P)
			end
			if x[i]>=0.5 && x[i]<1-P
				x[i+1]=(1-P-x[i])/(0.5-P)
			end
			if x[i]>=1-P && x[i]<1
				x[i+1]=(1-x[i])/P
			end    
			G[i]=x[i]*Value
		end

	elseif index == 7
		# Sine map
		for i=1:max_iter
			 x[i+1] = sin(pi*x[i])
			 G[i]=(x[i])*Value
		 end
	elseif index == 8
		 # Singer map 
		 u=1.07
		 for i=1:max_iter
			 x[i+1] = u*(7.86*x[i]-23.31*(x[i]^2)+28.75*(x[i]^3)-13.302875*(x[i]^4))
			 G[i]=(x[i])*Value
		 end
	elseif index == 9
		# Sinusoidal map
		 for i=1:max_iter
			 x[i+1] = 2.3*x[i]^2*sin(pi*x[i])
			 G[i]=(x[i])*Value
		 end
		 
	elseif index == 10
		 # Tent map
		 x[1]=0.6
		 for i=1:max_iter
			 if x[i]<0.7
				 x[i+1]=x[i]/0.7
			 end
			 if x[i]>=0.7
				 x[i+1]=(10/3)*(1-x[i])
			 end
			 G[i]=(x[i])*Value
		 end

	end
	return G[curr_iter]

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
		println("===========[CGSA results ]=============")
		println("| Generations = $max_it")
		println("| Evals       = ", max_it*N)
		@printf("| best f.     = %e\n", Fbest)
		println("=======================================")
	end

	return Lbest, Fbest

end