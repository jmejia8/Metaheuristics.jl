if VERSION >= v"0.7.0"
    
    using Pkg
    "Distributions" ∉ keys(Pkg.installed()) && Pkg.add("Distributions")

elseif VERSION < v"0.7.0" && Pkg.installed("Distributions") == nothing
    Pkg.add("Distributions")
end

using Distributions

if VERSION >= v"0.7.0"
	import Random: randperm
	import Printf.@printf
	import LinearAlgebra: norm, Diagonal, dot

	indmin = argmin
	indmax = argmax

	eye(n) = Diagonal(ones(n))

	repmat(X, n) = repeat(X, n)
end