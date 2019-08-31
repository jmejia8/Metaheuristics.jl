if VERSION >= v"0.7.0"
    
    using Pkg
    "Distributions" ∉ keys(Pkg.installed()) && Pkg.add("Distributions")

elseif VERSION < v"0.7.0" && Pkg.installed("Distributions") == nothing
    Pkg.add("Distributions")
end


if VERSION >= v"0.7.0"
	# import DelimitedFiles.writedlm
	import Random: randperm
	import Printf.@printf
	import LinearAlgebra: norm, Diagonal, dot
    import Statistics: var, mean, std

	indmin = argmin
	indmax = argmax

	eye(n) = Diagonal(ones(n))

	repmat(X, n) = repeat(X, n)
	# writecsv(f, A) = writedlm(f, A, ',')
end