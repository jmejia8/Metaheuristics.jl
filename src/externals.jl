if VERSION >= v"0.7.0"
	import Random: randperm, shuffle!, shuffle, seed!
	import Printf.@printf
	import LinearAlgebra
	import LinearAlgebra: norm, Diagonal, dot
	import Statistics: var, mean, std
	using UnicodePlots
	import Base.minimum
end

using Requires

function __init__()
	@require Optim ="429524aa-4258-5aef-a3af-852621145aeb" @eval import Optim
end

