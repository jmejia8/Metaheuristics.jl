"""
	spacing(A)
Computes the Schott spacing indicator. `spacing(A) == 0` means that vectors in `A` are
uniformly distributed.
"""
function spacing(A::Array{Array{Float64,1},1})
    n = length(A)

    d = zeros(n)
    for i = 1:n
    	@inbounds d[i] = minimum( [ sum( abs.(A[i] - A[j]) ) for j in deleteat!(collect(1:n), i) ] )
    end

    d_hat = mean(d)
    s = (1.0 / ( n - 1.0 )) * sum( (d_hat .- d) .^ 2 )
    return sqrt(s)
    
end


spacing(A::Array{T}) where T <: AbstractMultiObjectiveSolution = spacing( fval.(A) )
spacing(A::State{T}) where T <: AbstractMultiObjectiveSolution = spacing(A.population)
spacing(A::Matrix) = spacing([A[i,:] for i in 1:size(A,1)])
