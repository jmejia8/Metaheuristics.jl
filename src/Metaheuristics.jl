__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso, ecaConstrained, WOA, GOA, GSA

# ECA algorithm
include("eca.jl")
include("ecaConstrained.jl")

# Differential algorithm
include("diffEvolution.jl")

# PSO algorithm
include("pso.jl")

# WOA
include("WOA.jl")

# GOA
include("GOA.jl")

# GSA
include("GSA.jl")


end # module
