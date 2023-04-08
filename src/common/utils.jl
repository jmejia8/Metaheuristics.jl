function pairwise_distances(population::AbstractArray{T},args...;kargs...) where T <: AbstractSolution
    pairwise_distances(fvals(population), args...;kargs...)
end


function pairwise_distances(A::Matrix, dist = Euclidean(); diag_val=Inf)
    D = pairwise(dist, A, dims=1)
    D[Diagonal(ones(Bool, size(D,1)))] .= diag_val

    return D
end

_mat_to_bounds(bounds::AbstractMatrix) = BoxConstrainedSpace(lb=bounds[1,:], ub=bounds[2,:])
_mat_to_bounds(space::AbstractSearchSpace) =  space
