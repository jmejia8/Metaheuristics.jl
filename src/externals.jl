import Random: randperm, shuffle!, shuffle, seed!
import Printf.@printf
import LinearAlgebra
import LinearAlgebra: norm, Diagonal, dot
import Statistics: var, mean, std
import Distances: pairwise, evaluate, Euclidean, euclidean
using UnicodePlots
import Base.minimum

using Requires

function __init__()
	@require JMcDM ="358108f5-d052-4d0a-8344-d5384e00c0e5" include("DecisionMaking/JMcDM.jl")
end

