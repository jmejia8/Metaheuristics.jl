"""
RDEx: An effectiveness-driven hybrid Single-Objective optimizer

Julia implementation of the adaptive differential evolution algorithm.

"""
module RDExModule

# Standard library dependencies
using Random
using Statistics
using Distributions

# Submodules
include("RDEx/constants.jl")
include("RDEx/random.jl")
include("RDEx/selection.jl")
include("RDEx/types.jl")
include("RDEx/mutation.jl")
include("RDEx/adaptation.jl")
include("RDEx/metaheuristics.jl")


# Export public API
export RDEx



# Import Metaheuristics API (avoid conflict by importing functions directly)
import .MetaheuristicsAPI: RDEx 
# Note: RDExAlgorithm is the Metaheuristics-compatible algorithm constructor
# Users can do: using RDEx: RDExAlgorithm; algo = RDExAlgorithm()
# Or: import RDEx: RDExAlgorithm as RDEx; algo = RDEx()

end # module RDExModule
