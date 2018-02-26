include("tools.jl")
function IPD(N, D, XiMin, XiMax, γ)
	Population = zeros(N, D)
		
	# place probes probe line-by-probe line (i is dimension number)
	for i = 1:D 
		ΔXi = γ * (XiMax[i] - XiMin[i]) / (N - 1)
		Population[:, i] = XiMin[i] + ΔXi*linspace(0, N-1,N)
	end

	return Population

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

function CFO(fobj::Function,
                D::Int;
                N::Int  = 8D,
        max_evals::Int  = 10000D, 
                γ::Real = 1,
                α::Real = 2,
                β::Real = 2,
         FrepInit::Real = 0.5,
            ΔFrep::Real = 0.005,
          minFrep::Real = 0.05,
                limits = [-100.0, 100.0])

	XiMin, XiMax = limits
	XiMin *= ones(D)
	XiMax *= ones(D)

	max_iters = div(max_evals, N) + 1

	LastStep = max_iters

	# STEP (A1) -------- Compute Initial Probe Distribution (Step 0)-----------
	Population = initializePop(N, D, XiMin, XiMax)#IPD(N, D, XiMin, XiMax, γ)
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
			for i = 1:D
				for k = 1:N
					if k == p
						continue
					end

					sum_kp = sum((Population[k, :] - Population[p, :]) .^ 2)

					# to avoid zero denominator
					if sum_kp != 0 
						denom     = sqrt(sum_kp)
						Numerator = UnitStep(fitness[k] - fitness[p]) * (fitness[k] - fitness[p])
						accel[p, i] += (Population[k, i] - Population[p, i]) * (Numerator ^ α) / (denom ^ β)
					end
				end
			end
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
		Population += accel
		
		# Starting at Step #20 Shrink Decision Space Around Best Probe Every
		# 20th Step
		if mod(j, 20) == 0 && j >= 21
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

	return best, bestFitness
end
