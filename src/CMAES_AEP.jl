# Based on:
# Li, Zhenhua, and Qingfu Zhang. 
# "An efficient rank-1 update for Cholesky CMA-ES using auxiliary evolution path." 
# Evolutionary Computation (CEC), 2017 IEEE Congress on. IEEE, 2017.

using Distributions

struct IndCMA
	# step
	z::Vector
	y::Vector

	# position
	x::Vector

	# cost
	f::Real
end

function mvnorm(D)
	return rand(MvNormal(zeros(D), eye(D)))
end

function CMAES_AEP(fobj::Function,
                D::Int;
        max_evals::Int = 10000D, 
      showResults::Bool = true,
         saveLast::String = "",
  saveConvergence::String = "",
             limits = [-100.0, 100.0])

	# algorithm parameters
	λ = 4 + floor(Int, 3log(D))
	μ = div(λ, 2)
	w = (log(μ+1) - log.(1:μ)) ./ (μ*log(μ+1) - sum( log.(1:μ) ))
	μw = 1 / sum( w .^ 2 )
	cσ = √(μw) / ( √(D) + √(μw))
	dσ = 1 + 2max(0, √((μw - 1) / (D+1)) - 1) + cσ
	cc = 4 / (D + 4)
	c1 = 2 / ( D + √(2) )^2


	ENN= √(D)*(1-1/(4*D)+1/(21*D^2))

	# Limits
	VarMin, VarMax = limits


	σ = (VarMax - VarMin)/3

	# Auxiliary Evolution Path
	P = zeros(D)
	Pσ= zeros(D)
	V = zeros(D)

	# initialize M
	x = VarMin + (VarMax - VarMin) * rand(D)
	f = fobj(x)
	M = IndCMA(zeros(D), zeros(D), x, f)
	A = eye(D)

	# current number of evaluations
	nevals = 1

	# best solution
	bestSol = M

	# convergence values
	convergence = []
	if saveConvergence != ""
		push!(convergence, [nevals bestSol.f])
	end

	# stop condition
	stop = false

	t = 1
	while !stop

		# Generate Samples
		Population = Array{IndCMA}([])
		fVals = zeros(λ)
		for i=1:λ
			z = mvnorm(D)
			y = A*z
			x = M.x + σ * y
			f = fobj(x)

			nevals += 1

			sol = IndCMA(z, y, x, f)
			fVals[i] = f
			
			push!(Population, sol)

			# Update best solution
			if sol.f < bestSol.f
				bestSol = sol
			end
		end

		if saveConvergence != ""
			push!(convergence, [nevals bestSol.f])
		end

		stop = nevals >= max_evals

		if stop
			break
		end

		Population = Population[sortperm(fVals)]

		x   = zeros(D)
		y_w = zeros(D)
		z_w = zeros(D)

		for i = 1:μ
			x += w[i] * Population[i].x
			y_w += w[i] * Population[i].y
			z_w += w[i] * Population[i].z
		end

		M = IndCMA(zeros(D), zeros(D), x, 0)

		# update paths
		P = (1- cc)*P + √(cc*( 2 - cc )*μw) * y_w
		V = (1- cc)*V + √(cc*( 2 - cc )*μw) * z_w

		norm_v = dot(V, V)
		b = √(1-c1) / norm_v
		b*= √( 1 + norm_v * c1/(1-c1)) - 1

		A = √(1-c1) * A + b*P * V'

		Pσ = (1-cσ)*Pσ + √(cσ * (2-cσ) *μw )*z_w

		σ *= exp( (cσ/dσ)* (norm(Pσ) / ENN - 1) )
		t += 1
	end

	if saveConvergence != ""
		writecsv(saveConvergence, convergence)
	end

	if saveLast != ""
		writecsv(saveLast, M.x)        
	end

	if showResults
		println("===========[CMAES results]=============")
		println("| Generations = $t")
		println("| Evals       = ", nevals)
		@printf("| best f.     = %e\n", bestSol.f)
		println("=======================================")
	end

	return bestSol.x, bestSol.f

end