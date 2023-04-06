"""
    Grid(npartitions, dim)

Parameters to generate a grid with `npartitions` in a space with `dim` dimensions.

### Example 

```julia-repl
julia> sample(Grid(5,2))
25×2 Matrix{Float64}:
 0.0   0.0
 0.25  0.0
 0.5   0.0
 0.75  0.0
 ⋮     
 0.5   1.0
 0.75  1.0
 1.0   1.0

julia> sample(Grid(5,2), [-1 -1; 1 1.])
25×2 Matrix{Float64}:
 -1.0  -1.0
 -0.5  -1.0
  0.0  -1.0
  0.5  -1.0
  ⋮    
  0.0   1.0
  0.5   1.0
  1.0   1.0
```

Note that the sample is with size `npartitions^(dim)`.
"""
struct Grid <: AbstractInitializer
    npartitions::Int
    dim::Int
end

function sample(method::Grid, bounds = zeros(0,0))
    v = range(0, 1, length = method.npartitions)
    vals = Iterators.product((v for _ in 1:method.dim)...)
    
    # fill values in the grid
    _X = reshape([x for val in vals for x in val], method.dim, length(vals))
    X = Array(_X')

    isempty(bounds) && (return X)
    _scale_sample(X, bounds)
end


