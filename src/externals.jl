if VERSION >= v"0.7.0"
	# import DelimitedFiles.writedlm
	import Random: randperm, shuffle!
	import Printf.@printf
	import LinearAlgebra: norm, Diagonal, dot
    import Statistics: var, mean, std

	indmin = argmin
	indmax = argmax

	eye(n) = Diagonal(ones(n))

	repmat(X, n) = repeat(X, n)
	# writecsv(f, A) = writedlm(f, A, ',')
end
