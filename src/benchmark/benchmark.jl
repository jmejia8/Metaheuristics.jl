module Benchmark

include("box-constrained.jl")
include("constrained.jl")
include("multiobjective.jl")


function get_problem(problem)
    if problem == :sphere
        return sphere()
    elseif problem == :discus
        return discus()
    elseif problem == :rastrigin
        return rastrigin()
    end

    error("Argument should be: :sphere, :discus, :rastrigin")

end

end
