import Random: randperm, shuffle!, shuffle, seed!
import Printf.@printf
import LinearAlgebra
import LinearAlgebra: norm, Diagonal, dot
import Statistics: var, mean, std
import Distances: pairwise, evaluate, Euclidean, euclidean
import Base.minimum

using JMcDM
using Requires

using SearchSpaces
import SearchSpaces: AbstractSearchSpace, getdim

function __init__()
    @require UnicodePlots = "b8865327-cd53-5732-bb35-84acbb429228" begin
        include("common/show_plots.jl")
    end
end

