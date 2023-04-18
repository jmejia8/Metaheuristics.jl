"""
	covering(A, B)
Computes the covering indicator (percentage of vectors in B that are dominated by vectors in A)
from two sets with non-dominated solutions.

A and B with size (n, m) where n is number of samples and m is the vector dimension.

Note that `covering(A, B) == 1` means that all solutions in B are dominated by those
in A. Moreover, `covering(A, B) != covering(B, A)` in general.

If `A::State` and `B::State`, then computes `covering(A.population, B.population)` after
ignoring dominated solutions in each set.
"""
function covering(A::Array{Vector{T}}, B::Array{Vector{T}}) where T <: Real
    n_a = length(A)
    n_b = length(B)
    s = 0.0
    for i = 1:n_b
        for j = 1:n_a
			# `compare` returns 1 if argument 1 dominates argument 2
            if compare(A[j], B[i]) == 1
                s += 1.0
                break
            end    
        end
        
    end

    return s / n_b
    
end

covering(A::Vector{T}, B::Vector{T}) where T <: AbstractMultiObjectiveSolution = covering(fval.(A), fval.(B))

function covering(A::State{T}, B::State{T};
        verbose = true
    ) where T <: AbstractMultiObjectiveSolution

	A_non_dominated = get_non_dominated_solutions(A.population)
	B_non_dominated = get_non_dominated_solutions(B.population)
	
	length(A_non_dominated) != length(A.population) &&
	verbose && @warn "Some solutions in B dominate each other (solutions ignored)."


	length(B_non_dominated) != length(B.population) &&
	verbose && @warn "Some solutions in B dominate each other (solutions ignored)."

	covering(A_non_dominated, B_non_dominated)
end


function covering(A::Matrix, B::Matrix)
	A_arr = [ A[i,:] for i in 1:size(A,1) ]
	B_arr = [ B[i,:] for i in 1:size(B,1) ]
	covering(A_arr, B_arr)
end

