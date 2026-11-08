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

function sample(method::Grid, bounds = zeros(0, 0))
    n, d = method.npartitions, method.dim
    v = range(0, 1, length = n)
    m = 1
    for _ in 1:d
        m = Base.checked_mul(m, n) # number of points: n^d
    end

    # fill values in the grid (the first coordinate changes fastest)
    X = Matrix{Float64}(undef, m, d)
    @inbounds for row in 0:m - 1
        q = row
        for j in 1:d
            q, r = divrem(q, n)
            X[row + 1, j] = v[r + 1]
        end
    end

    isempty(bounds) && (return X)
    _scale_sample(X, bounds)
end
