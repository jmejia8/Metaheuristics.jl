_epsilon(A, B, j) = [maximum(A[i,:] ./ (B[j,:])) for i in 1:size(A, 1)]

"""
    epsilon_indicator(A, B)

Computes the ε-indicator for non-dominated sets `A` and `B`.
It is assumed that all values in `A` and `B` are positive. If negative,
the sets are translated to positive values.

### Interpretation

- `epsilon_indicator(A, PF)` is unary if `PF` is the Pareto-optimal front.
- `epsilon_indicator(A, B) == 1` none is better than the other.
- `epsilon_indicator(A, B) < 1` means that A is better than B.
- `epsilon_indicator(A, B) > 1` means that B is better than A.
- Values closer to 1 are preferable.

### Examples

```julia-repl
julia> A1 = [4 7;5 6;7 5; 8 4.0; 9 2];

julia> A2 = [4 7;5 6;7 5; 8 4.0];

julia> A3 = [6 8; 7 7;8 6; 9 5;10 4.0 ];

julia> PerformanceIndicators.epsilon_indicator(A1, A2)
1.0

julia> PerformanceIndicators.epsilon_indicator(A1, A3)
0.9

julia> f, bounds, pf = Metaheuristics.TestProblems.ZDT3();

julia> res = optimize(f, bounds, NSGA2());

julia> PerformanceIndicators.epsilon_indicator(res, pf)
1.00497701620997
```

"""
function epsilon_indicator(A::AbstractMatrix, B::AbstractMatrix)
    ma = minimum(A)
    mb = minimum(B)
    _c = min(ma, mb)
    
    mb == zero(typeof(mb)) && error("Second argument contains zeros.")
 
    if _c < zero(typeof(_c))
        A = A .- 2_c
        B = B .- 2_c
    end

    return maximum([minimum(_epsilon(A, B, j)) for j in 1:size(B,1)])
end

epsilon_indicator(A, B) = epsilon_indicator(fvals(A), fvals(B))

