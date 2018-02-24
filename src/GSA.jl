# GSA code v1.1.
# Coded by Jesús Mejúa, 2010. 
# Based on MATLAB code of Esmat Rashedi, 2010. 
# " E. Rashedi, H. Nezamabadi-pour and S. Saryazdi,
# “GSA: A Gravitational Search Algorithm”, Information sciences, vol. 179,
# no. 13, pp. 2232-2248, 2009."

function initializationGSA(D,N,up,down)
	if size(up,2) == 1
		X = rand(N,D).*(up-down)+down
	end
	if size(up,2)>1
		for i = 1:D
			high = up[i]
			low  = down[i]
			X[:,i] = rand(N,1).*(high-low)+low
		end
	end

	return X
end

function evaluateF(X, fobj, N)

	f = zeros(N)
	for i = 1:N
		f[i] = fobj(X[i,:])
	end

	return f
end

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
		
		M = (fit - worst) ./ (best - worst) #eq.15,

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
					E[i, k] = E[i, k]+rand() * (M[j]) * ((X[j,k]-X[i, k]) / (R^Rpower))
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

function getBestGSA(fitness, searchType)
	if searchType == :minimize
		best_X = indmin(fitness) # minimization.
		best = fitness[best_X] 
	else
		best_X = indmax(fitness) # maximization.
		best = fitness[best_X] 
	end

	return best_X, best
end

# Gravitational Search Algorithm.
function GSA(fobj::Function,
	            D::Int;
		        N::Int    = 30,
	 ElitistCheck::Int    = 1,
	   searchType::Symbol = :minimize,
		max_evals::Int    = 10000D,
		   Rpower::Int    = 1,
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

	max_it = div(max_evals, N) + 1

	#random initialization for agents.
	X = initializationGSA(D,N,up,low)

	#Evaluation of agents. 
	fitness = evaluateF(X, fobj, N)

	V = zeros(N,D)

	# Current best
	best_X, best = getBestGSA(fitness, searchType)

	Fbest = best
	Lbest = X[best_X,:]

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
		# X = space_bound(X,up,low)
		
		#Evaluation of agents. 
		fitness = evaluateF(X, fobj, N)
		best_X, best = getBestGSA(fitness, searchType)


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

	end #iteration

	return Lbest, Fbest

end