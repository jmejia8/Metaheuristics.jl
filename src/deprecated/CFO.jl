function IPD(N, D, XiMin, XiMax, γ)
	if D > 1
		NumProbesPerDimension = div(N, D); # even #
	else
		NumProbesPerDimension = N;
	end

	R = repmat(XiMin' + γ .* (XiMax' - XiMin'), N);

	for i = 1: D # place probes probe line-by-probe line (i is dimension number)
		DeltaXi = (XiMax[i] - XiMin[i]) / (NumProbesPerDimension - 1);
		for k = 1:NumProbesPerDimension
			p = k + NumProbesPerDimension * (i - 1); # probe #
			R[p, i] = XiMin[i] + (k - 1) * DeltaXi;
		end
	end

	return R

end

function RetrieveErrantProbes(Population, XiMin, XiMax, Frep)
	N = size(Population, 1)
	repXiMin = repmat(XiMin', N)
	repXiMax = repmat(XiMax', N)

	Rcomp = Population .< repXiMin
	Population = .!Rcomp .* Population + Rcomp .* max.(repXiMin + Frep .* (Population - repXiMin), repXiMin)

	Rcomp = Population .> repXiMax
	Population = .!Rcomp .* Population + Rcomp .* min.(repXiMax - Frep .* (repXiMax - Population), repXiMax)

	return Population
end

function UnitStep(vars)
	return vars >= 0

end

function CFO(fobj_::Function,
                D::Int;
        max_evals::Int  = 10000D, 
                γ::Real = 0.5,
                G::Real = 1,
                α::Real = 1,
                β::Real = 2,
         FrepInit::Real = 0.5,
            ΔFrep::Real = 0.1,
          minFrep::Real = 0.05,
       searchType::Symbol=:minimize,
                limits = [-100.0, 100.0])

	if searchType == :minimize
		fobj(x) = -fobj_(x)
	else
		fobj = fobj_
	end

	XiMin, XiMax = limits
	XiMin *= ones(D)
	XiMax *= ones(D)

	if D == 30
		ρ = 2
	elseif D <= 30
		ρ = 12
	else
		ρ = 8
	end

	N = ρ*D

	max_iters = div(max_evals, N) + 1

	LastStep = max_iters

	# STEP (A1) -------- Compute Initial Probe Distribution (Step 0)-----------
	Population = IPD(N, D, XiMin, XiMax, γ)
	# return Population

	# STEP (A2) ---------- Compute Initial Fitness Matrix (Step 0) ------------
	fitness = zeros(N)
	fitness = evaluatePop(Population, fobj, N)
	neval = N

	# STEP (A3) ------- Set Initial Probe Accelerations to ZERO (Step 0)-------
	accel = zeros(N, D)

	# STEP (A4) --------------------- Initialize Frep -------------------------
	Frep = FrepInit


	# ======= LOOP ON TIME STEPS STARTING AT STEP #1 (#2 in this code) ========

	best_i, bestFitness = getBest(fitness, :maximize)
	best = Population[best_i, :]

	for j = 1:max_iters
		
		# Compute Accelerations Based on Current Probe Distribution & Fitnesses
		# ---------------------------------------------------------------------
		for p = 1:N
			for k = 1:N
				if k == p
					continue
				end

				sum_kp = sum((Population[k, :] - Population[p, :]) .^ 2)

				# to avoid zero denominator
				if sum_kp != 0 
					denom     = sqrt(sum_kp)
					Numerator = UnitStep(fitness[k] - fitness[p]) * (fitness[k] - fitness[p])^α
					accel[p,:] += (Population[k,:] - Population[p, :]) * (Numerator) / (denom ^ β)
				end
			end

			accel[p,:] *= G
		end
		
		# Get Best Fitness & Corresponding Probe # and Time Step -------
		best_i, currentBestF = getBest(fitness, :maximize)
		
		if currentBestF >= bestFitness
			bestFitness = currentBestF
			best = Population[best_i, :]
		end
		
		# ------------------------ Increment Frep -----------------------------
		Frep = Frep + ΔFrep
		if Frep >= (1 + eps()) # keep Frep in range [0.05,1]
			Frep = minFrep
		end
		
		# STEP (B) ---- Compute Probe Position Vectors for this Time Step -----
		Population += 0.5accel
		
		# Starting at Step #20 Shrink Decision Space Around Best Probe Every
		# 20th Step
		if mod(j, 20) == 0
			# shrink DS by 0.5
			XiMin = XiMin + (best - XiMin) ./ 2
			XiMax = XiMax - (XiMax - best) ./ 2			
		end

		
		# Retrieve Errant Probes --------------------
		Population = RetrieveErrantProbes(Population, XiMin, XiMax, Frep)
		

		# Compute Fitness for Current Probe Distribution --
		fitness = evaluatePop(Population, fobj, N)
		neval += N
		
	end

	if searchType == :minimize
		bestFitness *= -1
	end
	return best, bestFitness
end
