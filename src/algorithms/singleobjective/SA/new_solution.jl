
function newSol(y, μ)
	# This function is used to generate new point according to lower and upper
	# and a random factor proportional to current point.
	return(((1.0 .+ μ) .^ abs.(y) .- 1) / μ) .* sign.(y)
end


