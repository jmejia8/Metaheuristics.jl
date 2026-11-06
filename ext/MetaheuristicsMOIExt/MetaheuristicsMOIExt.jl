module MetaheuristicsMOIExt

using Metaheuristics
using SearchSpaces
import MathOptInterface as MOI

function __init__()
    @static if isdefined(Base, :get_extension)
        setglobal!(Metaheuristics, :Optimizer, Optimizer)
    end
    println("Extension was loaded!")
end

include("constants.jl")
include("optimizer.jl")
include("utils.jl")
include("legacy.jl")
include("initialization.jl")
include("multiobjective.jl")
include("optimize.jl")
include("results.jl")

end # module
