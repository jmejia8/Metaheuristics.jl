"""
SA:  Minimizes a function with the method of Simulated Annealing
(Kirkpatrick et al., 1983)

"""

function newSol(y, μ)
	# This function is used to generate new point according to lower and upper
	# and a random factor proportional to current point.
	return(((1+ μ).^abs.(y)-1)/ μ).*sign.(y)
end


function SA( fobj::Function,
                D::Int;
               x0::Vector = [],
                N::Int = 500,
        max_evals::Int = 10000D,
           TolFun::Real= 1e-4,
      showResults::Bool= true,
               limits  = [-100., 100.])

	l, u = limits

	if length(x0) != D
	    x0 = l + (u - l) .* rand(D)
	end

	# the current point and fx=f(x)
	x = x0
	fx= fobj(x)
	f0= fx


	nevals = 1
	stop = false

	# Main loop simulates de annealing from a high temperature to zero in max_evals.
	while !stop

		# T is the inverse of temperature.
		T = nevals / max_evals 
		μ = 10.0 ^( 100T )    
		
		# For each temperature we take 500 test points to simulate reach termal
		# equilibrium.
		for i = 1:N        
			# We generate new test point using newSol function      
			dx = newSol(2rand(size(x)) - 1, μ) .* (u-l)

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
	end
	
	return x0, f0
end
