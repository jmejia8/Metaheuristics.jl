include("gen_initial_state.jl")

"""
    set_user_solutions!(optimizer, x, fx;verbose=true)

Provide initial solutions to the `optimizer`.
- `x` can be a `Vector` and `fx` a function or `fx = f(x)`
- `x` can be a matrix containing solutions in rows.

### Example

```julia-repl
julia> f(x) = abs(x[1]) + x[2]  + x[3]^2 # objective function
f (generic function with 1 method)

julia> algo  = ECA(N = 61); # optimizer

julia> # one solution can be provided
       x0 = [0.5, 0.5, 0.5];

julia> set_user_solutions!(algo, x0, f);

julia> # providing multiple solutions
       X0 = rand(30, 3); # 30 solutions with dim 3

julia> set_user_solutions!(algo, X0, f);

julia> optimize(f, [0 0 0; 1 1 1.0], algo)
Optimization Result
===================
  Iteration:       235
  Minimum:         1.69759e-46
  Minimizer:       [6.87417e-47, 7.75848e-47, 4.84071e-24]
  Function calls:  14305
  Total time:      0.1050 s
  Stop reason:     Due to Convergence Termination criterion.
```

"""
function set_user_solutions!(algo::AbstractAlgorithm, solution::AbstractSolution;kargs...)
    status = algo.status
    if !isnothing(status.best_sol) && !isempty(status.population)
        push!(status.population, solution)
    else
        algo.status = State(solution, [solution])
    end
    algo
end


function set_user_solutions!(algo::AbstractAlgorithm, x::AbstractVector, fx;kargs...)
    set_user_solutions!(algo, create_child(x, fx); kargs...)
end

function set_user_solutions!(algo::AbstractAlgorithm, x::AbstractVector, f::Function;kargs...)
    set_user_solutions!(algo, x, f(x); kargs...)
end


function set_user_solutions!(algo::AbstractAlgorithm, X::AbstractMatrix, fX::AbstractVector;verbose=true)
    n = size(X, 1)
    m = length(fX)

    if n != m
        verbose && @warn "$(n) decision vectors provided but $(m) objective values."
        n = min(m, n)
        verbose && println("Taking ", n, " as the number of initial solutions.")
    end

    # nothing to do due to it is necessary the objective value
    n == 0 && (return algo)

    # TODO: this part can be parallelized
    for i in 1:n
        set_user_solutions!(algo, X[i,:], fX[i])
    end

    # TODO check population size provided in algo.parameters.N

    algo
end

function set_user_solutions!(algo::AbstractAlgorithm, X::AbstractMatrix, f::Function;kargs...)
    set_user_solutions!(algo, X, [f(X[i,:]) for i in 1:size(X,1)];kargs...)
end

