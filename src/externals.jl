import Random
import Random: randperm, shuffle!, shuffle, seed!
import Printf: @printf, @sprintf
import LinearAlgebra
import LinearAlgebra: norm, Diagonal, dot
import Statistics: var, mean, std
import Distances: pairwise, evaluate, Euclidean, euclidean
import Base.minimum

using JMcDM
using Requires
using Reexport

@reexport using SearchSpaces
import SearchSpaces: AbstractSearchSpace,getdim
import SearchSpaces: AtomicSearchSpace, AbstractSampler, Sampler

function __init__()
    @require UnicodePlots = "b8865327-cd53-5732-bb35-84acbb429228" begin
        include("common/show_plots.jl")
    end
end

