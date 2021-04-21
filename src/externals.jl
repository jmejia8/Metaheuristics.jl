if VERSION >= v"0.7.0"
	import Random: randperm, shuffle!, shuffle, seed!
	import Printf.@printf
	import LinearAlgebra: norm, Diagonal, dot
	import Statistics: var, mean, std
	using UnicodePlots
	import Base.minimum
end
