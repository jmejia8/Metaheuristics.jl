"""
	deltap(front, true_pareto_front; p = 1)
	Δₚ(front, true_pareto_front; p = 1) 

Returns the averaged Hausdorff distance indicator aka Δₚ (Delta p).

"`Δₚ`" can be typed as `\\Delta<tab>\\_p<tab>`.

### Parameters

`front` and `true_pareto_front` can be:

- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `Array{xFgh_indiv}` (usually `State.population`)

"""
function deltap(front::Array{Vector{T}}, true_pareto_front::Array{Vector{T}}; p = 1) where T <: Real
	distances = zeros(length(true_pareto_front), length(front))

	for i in eachindex(true_pareto_front)
		for j in eachindex(front)
			distances[i,j] = norm(true_pareto_front[i] - front[j])
		end
	end

	gd  = mean(minimum(distances, dims = 1) .^ p)
	igd = mean(minimum(distances, dims = 2) .^ p)
	return max(gd^(1/p), igd^(1/p))
end


deltap(front::Vector{T}, true_pareto_front::Vector{T}; kargs...) where T <: AbstractMultiObjectiveSolution = deltap(fval.(front), fval.(true_pareto_front); kargs...)

function deltap(front::AbstractMatrix, true_pareto_front::AbstractMatrix; kargs...)
	front_ = [ front[i,:] for i in 1:size(front,1) ]
	true_pareto_front_ = [ true_pareto_front[i,:] for i in 1:size(true_pareto_front,1) ]
	deltap(front_, true_pareto_front_; kargs...)
end

const Δₚ = deltap
