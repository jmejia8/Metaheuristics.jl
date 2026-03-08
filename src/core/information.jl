struct Information
    f_optimum::Float64
    x_optimum::Array{Float64}
end

"""
    Information Structure

`Information` can be used to store the true optimum in order to stop a metaheuristic early.

Properties:

- `f_optimum` known minimum.
- `x_optimum` known minimizer.

If `Options` is provided, then `optimize` will stop when `|f(x) - f(x_optimum)| < Options.f_tol`
or `‖ x - x_optimum ‖ < Options.x_tol` (euclidean distance).

# Example

If you want an approximation to the minimum with accuracy of `1e-3` (|f(x) - f(x*)| < 1e-3),
then you may use `Information`.

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> information = Information(f_optimum = 0.0)
Information(0.0, Float64[])

julia> options = Options(f_tol = 1e-3, seed=1)
Options(0.0, 0.001, 0.0, 0.0, 1000.0, 0.0, 0.0, 1, false, true, false, :minimize)

julia> state = optimize(f, bounds, ECA(information=information, options=options))
Optimization Result
===================
  Iteration:       25
  Minimum:         0.000450013
  Minimizer:       [0.0142217, 0.010249, -0.0162118]
  Function calls:  525
  Total time:      0.0263 s
  Stop reason:     The desired accuracy was obtained.
```
"""
function Information(;#
    f_optimum = NaN,
    x_optimum::Array{Float64} = Float64[],
)

    Information(Float64(f_optimum), x_optimum)

end

