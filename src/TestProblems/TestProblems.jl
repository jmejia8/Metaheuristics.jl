module TestProblems

import ..generateChild

include("box-constrained.jl")
include("constrained.jl")
include("multiobjective.jl")


function get_problem(problem)
    expr = Meta.parse( string(problem) * "()")
    return eval(expr)
    try
        return eval(expr)
    catch
    end
    error("Argument should be: :sphere, :discus, :rastrigin")

end

end
