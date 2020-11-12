# GSA code v0.1.
# Coded by Jesús Mejía. 
# Based on MATLAB code of Esmat Rashedi, 2010. 
# " E. Rashedi, H. Nezamabadi-pour and S. Saryazdi,
# “GSA: A Gravitational Search Algorithm”, Information sciences, vol. 179,
# no. 13, pp. 2232-2248, 2009."

function massCalculation(fit,searchType)
	####here, make your own function of 'mass calculation'

	Fmax = maximum(fit)
	Fmin = minimum(fit)
	Fmean= mean(fit) 
	N = length(fit)

	if Fmax == Fmin
		M = ones(N)
	else
		
		if searchType == :minimize #for minimization
			best = Fmin
			worst= Fmax #eq.17-18.
		else #for maximization
			best = Fmax
			worst= Fmin #eq.19-20.
		end
		
		M = (fit .- worst) ./ (best .- worst) #eq.15,

	end

	M = M ./ sum(M) #eq. 16.
	
	return M

end

function Gconstant(iteration,max_it)
	# here, make your own function of 'G'
	α = 20
	G0= 100

	G = G0*exp(-α * iteration / max_it) #eq. 28.

	return G
end

#This function calculates the accelaration of each agent in gravitational field. eq.7-10,21.
function Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)
	N, D   = size(X)
	final_per= 2 #In the last iteration, only 2 percent of agents apply force to the others.

	####total force calculation
	if ElitistCheck == 1
		kbest = final_per+(1-iteration/max_it)*(100-final_per) #kbest in eq. 21.
		kbest = round(Int, N*kbest/100)
	else
		kbest = N #eq.9.
	end
	
	ds = sortperm(M, rev=true)

	E = zeros(N, D)
	for i=1:N
		for ii=1:kbest
			j = ds[ii]
			if j != i
				R = norm(X[i,:]-X[j,:], Rnorm) #Euclidian distanse.
				for k=1:D 
                    E[i, k]= E[i,k]+rand() * (M[j]) * ((X[j,k]-X[i,k])/(R^Rpower + eps()))
					#note that Mp(i)/Mi(i)=1
				end
			end
		end
	end

	##acceleration
	a = E.*G #note that Mp(i)/Mi(i)=1

	return a
end

# This function updates the velocity and position of agents.
function move(X,a,V)
	# movement.
	N, D = size(X)
	V = rand(N, D).*V + a # eq. 11.
	X = X + V # eq. 12.

	return X, V
end

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
