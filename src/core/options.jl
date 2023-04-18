default_rng_mh(seed=1223334444) = Random.default_rng(Int(seed))

#####################################################
#
#         STRUCTURES FOR THE OPTIONS
#
#####################################################
mutable struct Options
    x_tol::Float64
    f_tol::Float64
    g_tol::Float64
    h_tol::Float64
    f_tol_rel::Float64
    f_calls_limit::Float64
    time_limit::Float64

    iterations::Int
    store_convergence::Bool
    debug::Bool
    seed::UInt
    parallel_evaluation::Bool
    verbose::Bool
    rng
end

"""
    Options(;
        x_tol::Real = 1e-8,
        f_tol::Real = 1e-12,
        f_tol_rel::Real = eps(),
        f_tol_abs::Real = 0.0,
        g_tol::Real = 0.0,
        h_tol::Real = 0.0,
        f_calls_limit::Real = 0,
        time_limit::Real = Inf,
        iterations::Int = 1,
        store_convergence::Bool = false,
        debug::Bool = false,
        seed = rand(UInt),
        parallel_evaluation = false,
        verbose = false,
    )


`Options` stores common settings for metaheuristics such as the maximum number of iterations
debug options, maximum number of function evaluations, etc.

Main properties:

- `x_tol` tolerance to the true minimizer if specified in `Information`.
- `f_tol` tolerance to the true minimum if specified in `Information`.
- `f_tol_rel` relative tolerance.
- `f_calls_limit` is the maximum number of function evaluations limit.
- `time_limit` is the maximum time that `optimize` can spend in seconds.
- `iterations` is the maximum number of allowed iterations.
- `store_convergence` if `true`, then push the current `State` in `State.convergence` at each generation/iteration
- `debug` if `true`, then `optimize` function reports the current `State` (and interest information) for each iterations.
- `seed` non-negative integer for the random generator seed.
- `parallel_evaluation` enables batch evaluations.
- `verbose` show simplified results each iteration.
- `rng` user-defined Random Number Generator.

# Example

```julia-repl
julia> options = Options(f_calls_limit = 1000, debug=false, seed=1);

julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> state = optimize(f, bounds, ECA(options=options));
```

"""
function Options(;
        x_tol = 1e-8,
        f_tol = 1e-12,
        g_tol = 0.0,
        h_tol = 0.0,
        f_tol_rel = eps(),
        f_tol_abs = 0.0,
        f_calls_limit = 0.0,
        g_calls_limit = 0.0,
        h_calls_limit = 0.0,
        iterations::Int = 0,
        time_limit::Float64 = Inf,
        store_convergence = false,
        show_results = true,
        debug = false,
        search_type = :minimize,
        parallel_evaluation = false,
        verbose = false,
        seed = rand(UInt32),
        rng  = default_rng_mh(seed),
    )


    Options(
            promote(x_tol, f_tol, g_tol, h_tol, f_tol_rel)...,
            promote(f_calls_limit, time_limit)...,
            promote(iterations)...,
            # Results options
            promote(store_convergence, debug)...,
            UInt(seed),
            Bool(parallel_evaluation),
            Bool(verbose),
            rng
           )

end


function Base.show(io::IO, options::Options)
    _print_title(io, "Options")

    txt = String[]
    ln = Int[]
    for field in fieldnames(Options)
        v = getfield(options, field)
        fl = string(field) .* ":"
        push!(txt, @sprintf("  %-16s %s\n", fl, string(v)))
        push!(ln, length(fl))
    end
    println(io, join(txt[sortperm(ln)]))
end
