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

@reexport import SearchSpaces
@reexport import SearchSpaces: BoxConstrainedSpace, PermutationSpace, BitArraySpace, MixedSpace
@reexport import SearchSpaces: cardinality
using SearchSpaces
import SearchSpaces: AbstractSearchSpace,getdim
import SearchSpaces: AtomicSearchSpace, AbstractSampler, Sampler

"""
    boxconstraints(lb, ub, [rigid])
    BoxConstrainedSpace(ub, lb, [rigid])

Define a box-constrained search space (alias for BoxConstrainedSpace).

See [`SearchSpaces.BoxConstrainedSpace`](@ref) for more details.

# Example

```julia
f(x) = (x[1] - 100)^2 + sum(abs.(x[2:end]))
bounds = boxconstraints(lb = zeros(5), ub = ones(5), rigid = false)
optimize(f, bounds, ECA)
```
"""
const boxconstraints = BoxConstrainedSpace

function __init__()
    @require UnicodePlots = "b8865327-cd53-5732-bb35-84acbb429228" begin
        include("common/show_plots.jl")
    end
end

