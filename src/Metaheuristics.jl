__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso

# ECA algorithm
include("eca.jl")

# Differential algorithm
include("diffEvolution.jl")

# PSO algorithm
include("pso.jl")

end # module
