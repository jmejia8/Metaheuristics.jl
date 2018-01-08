__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso, ecaConstrained

# ECA algorithm
include("eca.jl")
include("ecaConstrained.jl")

# Differential algorithm
include("diffEvolution.jl")

# PSO algorithm
include("pso.jl")

end # module
