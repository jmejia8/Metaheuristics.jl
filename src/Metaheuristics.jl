__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso, ecaConstrained, WOA, GOA, GSA, SA

include("tools.jl")

# ECA algorithm
include("eca.jl")
include("ecaConstrained.jl")

# Differential algorithm
include("diffEvolution.jl")

# PSO algorithm
include("pso.jl")

# The whale optimization algorithm
# S Mirjalili, A Lewis - Advances in Engineering Software, 2016
include("WOA.jl")

# Grasshopper optimisation algorithm: Theory and application
# S Saremi, S Mirjalili, A Lewis - Advances in Engineering Software, 2017
include("GOA.jl")

# GSA: a gravitational search algorithm
# E Rashedi, H Nezamabadi-Pour, S Saryazdi - Information sciences, 2009
include("GSA.jl")

include("SA.jl")


end # module
