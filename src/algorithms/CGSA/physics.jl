
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
