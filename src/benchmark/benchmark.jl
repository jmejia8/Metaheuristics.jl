module Benchmark

include("box-constrained.jl")
include("constrained.jl")
include("multiobjective.jl")


function get_problem(problem)
    expr = Meta.parse( string(problem) * "()")
    try
        return eval(expr)
    catch
    end
    error("Argument should be: :sphere, :discus, :rastrigin")

end

end
