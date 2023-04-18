function pairwise_distances(population::AbstractArray{T},args...;kargs...) where T <: AbstractSolution
    pairwise_distances(fvals(population), args...;kargs...)
end


function pairwise_distances(A::Matrix, dist = Euclidean(); diag_val=Inf)
    D = pairwise(dist, A, dims=1)
    D[Diagonal(ones(Bool, size(D,1)))] .= diag_val

    return D
end

function _mat_to_bounds(bounds::AbstractMatrix)
    # validate bounds in cols
    if size(bounds, 1) > 2 && size(bounds, 2) == 2
        bounds = bounds'
    elseif size(bounds, 1) != 2
        error("Provide valid bounds. Suggestion: set `bounds = BoxConstrainedSpace(lb, ub)`.")
    end
    
    lb = bounds[1,:]
    ub = bounds[2,:]
    BoxConstrainedSpace(lb, ub)
end
_mat_to_bounds(space::AbstractSearchSpace) =  space
function _mat_to_bounds(b::Tuple{AbstractVector, AbstractVector})
    @assert length(b) == 2 "Provide valid lower and upper bounds."
    BoxConstrainedSpace(first(b), last(b))
end

