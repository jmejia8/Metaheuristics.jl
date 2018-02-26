function IPD(N, D, max_iters, XiMin, XiMax, γ, InitialProbesPlaceMode)

	R = zeros(N, D, max_iters)

	
	if InitialProbesPlaceMode == "Uniform On-Axis"
		if D > 1
			NumProbesPerDimension = div(N, D) # even #
		else
			NumProbesPerDimension = N
		end
		
		R[:, :, 1] = repmat(XiMin + γ .* (XiMax - XiMin), N)
		
		# place probes probe line-by-probe line (i is dimension number)
		for i = 1: D 
			DeltaXi = (XiMax[i] - XiMin[i]) / (NumProbesPerDimension - 1)
			for k = 1:NumProbesPerDimension
				p = k + NumProbesPerDimension * (i - 1) # probe #
				R[p, i, 1] = XiMin[i] + (k - 1) * DeltaXi
			end
		end
		
	elseif InitialProbesPlaceMode == "2D Grid"
		NumProbesPerDimension = sqrt(N)
		NumX1points = NumProbesPerDimension
		NumX2points = NumX1points
		DelX1 = (XiMax[1] - XiMin[1]) / (NumX1points - 1)
		DelX2 = (XiMax[2] - XiMin[2]) / (NumX2points - 1)
		
		for x1pointNum = 1:NumX1points
			for x2pointNum = 1:NumX2points
				p = NumX1points * (x1pointNum - 1) + x2pointNum  # probe #
				R[p, 1, 1] = XiMin[1] + DelX1 * (x1pointNum - 1) # x1 coord
				R[p, 2, 1] = XiMin[2] + DelX2 * (x2pointNum - 1) # x2 coord
			end
		end
		
	elseif InitialProbesPlaceMode == "Uniform On-Diagonal"
		DeltaX = (XiMax - XiMin) ./ (N - 1)
		R[:, :, 1] = repmat(XiMin, N, 1) + ((1:N) - 1)' * DeltaX
		
	elseif InitialProbesPlaceMode == "Random"
		R[:, :, 1] = repmat(XiMin, N, 1) + repmat((XiMax - XiMin), N, 1) .* rand(N, D)
	else
		error("Invalid InitialProbesPlaceMode.")
	end

	return R

end

function RetrieveErrantProbes(R, j, XiMin, XiMax, Frep)
	N = size(R, 1)
	repXiMin = repmat(XiMin', N)
	repXiMax = repmat(XiMax', N)

	Rcomp = R[:, :, j] .< repXiMin
	R[:, :, j] = .!Rcomp .* R[:, :, j] + Rcomp .* max.(repXiMin + Frep .* (R[:, :, j - 1] - repXiMin), repXiMin)

	Rcomp = R[:, :, j] .> repXiMax
	R[:, :, j] = .!Rcomp .* R[:, :, j] + Rcomp .* min.(repXiMax - Frep .* (repXiMax - R[:, :, j - 1]), repXiMax)

	return R
end

function UnitStep(vars)
	return vars >= 0

end

function HasFITNESSsaturated(Nsteps, j, M)
	out = false

	# execute at least 10 steps after averaging interval before performing this check
	if j < (Nsteps + 11) 
		return out
	end

	# tolerance for FITNESS saturation
	FitnessSatTOL = 0.000001 
	MeanFitness = mean(maximum(M[:, (j - Nsteps + 1): j], 1))
	BestFitness = maximum(M[:, j])

	# saturation if (avg value - last value) are within TOL
	if abs(MeanFitness - BestFitness) <= FitnessSatTOL 
		out = true
	end
	return out

end



function CFO(fobj::Function,
				D::Int;
				N::Int  = 8D,
        max_evals::Int  = 10000D, 
				γ::Real = 2,
				α::Real = 2,
				β::Real = 2,
		 FrepInit::Real = 0.5,
		    ΔFrep::Real = 0.005,
		  minFrep::Real = 0.05,
     InitialProbesPlaceMode = "Uniform On-Axis",  
				limits = [-100, 100])

	XiMin, XiMax = limits
	XiMin *= ones(D)
	XiMax *= ones(D)

	max_iters = div(max_evals, N) + 1

	LastStep = max_iters

	# STEP (A1) -------- Compute Initial Probe Distribution (Step 0)-----------
	R = IPD(N, D, max_iters, XiMin, XiMax, γ, InitialProbesPlaceMode)

	# STEP (A2) ---------- Compute Initial Fitness Matrix (Step 0) ------------
	M = zeros(N, max_iters)
	M[:, 1] = fobj(R[:, :, 1])
	Neval = N

	# STEP (A3) ------- Set Initial Probe Accelerations to ZERO (Step 0)-------
	A = zeros(N, D, max_iters)

	# STEP (A4) --------------------- Initialize Frep -------------------------
	Frep = FrepInit


	# ======= LOOP ON TIME STEPS STARTING AT STEP #1 (#2 in this code) ========
	BestFitness     = M[1, 1]
	BestProbeNumber = 1
	BestTimeStep    = 0

	for j = 2:max_iters
		
		# STEP (B) ---- Compute Probe Position Vectors for this Time Step -----
		R[:, :, j] = R[:, :, j - 1] + A[:, :, j - 1]
		
		# STEP (C) ---------------- Retrieve Errant Probes --------------------
		R = RetrieveErrantProbes(R, j, XiMin, XiMax, Frep)
		
		# STEP (D) --- Compute Fitness Matrix for Current Probe Distribution --
		M[:, j] = fobj(R[:, :, j])
		Neval = Neval + N
		
		# STEP (E) ------------------------------------------------------------
		# Compute Accelerations Based on Current Probe Distribution & Fitnesses
		# ---------------------------------------------------------------------
		for p = 1:N
			for i = 1:D
				for k = 1:N
					if k == p
						continue
					end

					SumSQ = sum((R[k, :, j] - R[p, :, j]) .^ 2)
					# to avoid zero denominator
					if SumSQ != 0 
						Denom = sqrt(SumSQ)
						Numerator = UnitStep(M[k, j] - M[p, j]) * (M[k, j] - M[p, j])
						A[p, i, j] = A[p, i, j] + (R[k, i, j] - R[p, i, j]) * (Numerator ^ α) / (Denom ^ β)
					end
				end
			end
		end
		
		# ------ Get Best Fitness & Corresponding Probe # and Time Step -------
		tmpBestProbeNumber = indmax(M[:, j])
		tmpBestFitness     = M[tmpBestProbeNumber, j] 
		
		if tmpBestFitness >= BestFitness
			BestFitness = tmpBestFitness
			BestProbeNumber = tmpBestProbeNumber
			BestTimeStep = j - 1
		end
		
		# ------------------------ Increment Frep -----------------------------
		Frep = Frep + ΔFrep
		if Frep >= (1 + eps()) # keep Frep in range [0.05,1]
			Frep = minFrep
		end
		
		# Starting at Step #20 Shrink Decision Space Around Best Probe Every
		# 20th Step
		if mod(j - 1, 20) <= eps() && j >= 21
			# shrink DS by 0.5
			XiMin = XiMin + (R[BestProbeNumber, :, BestTimeStep + 1] - XiMin) ./ 2
			XiMax = XiMax - (XiMax - R[BestProbeNumber, :, BestTimeStep + 1]) ./ 2
			
			R = RetrieveErrantProbes(R, j, XiMin, XiMax, Frep)
			# TO RETRIEVE PROBES LYING OUTSIDE SHRUNKEN DS
		end
		
		# STEP (F) ------------- Check for Early Run Termination --------------
		if HasFITNESSsaturated(25, j, M)
			LastStep = j - 1
			break
		end
	end

	return BestFitness#, BestProbeNumber, BestTimeStep, LastStep, R, M, Neval
end
