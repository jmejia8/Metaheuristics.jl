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
    f_calls_limit::Float64
    g_calls_limit::Float64
    h_calls_limit::Float64
    time_limit::Float64

    iterations::Int
    store_convergence::Bool
    show_results::Bool
    debug::Bool
    search_type::Symbol
    seed::UInt
    parallel_evaluation::Bool
    rng
end

"""
    Options(;
        x_tol::Real = 0.0,
        f_tol::Real = 0.0,
        g_tol::Real = 0.0,
        h_tol::Real = 0.0,
        f_calls_limit::Real = 0,
        g_calls_limit::Real = 0,
        h_calls_limit::Real = 0,
        time_limit::Real = Inf,
        iterations::Int = 1,
        store_convergence::Bool = false,
        debug::Bool = false,
        seed = rand(UInt),
        parallel_evaluation = false,
    )


`Options` stores common settings for metaheuristics such as the maximum number of iterations
debug options, maximum number of function evaluations, etc.

Main properties:

- `x_tol` tolerance to the true minimizer if specified in `Information`.
- `f_tol` tolerance to the true minimum if specified in `Information`.
- `f_calls_limit` is the maximum number of function evaluations limit.
- `time_limit` is the maximum time that `optimize` can spend in seconds.
- `iterations` is the maximum number iterationn permited.
- `store_convergence` if `true`, then push the current `State` in `State.convergence` at each generation/iteration
- `debug` if `true`, then `optimize` function reports the current `State` (and interest information) for each iterations.
- `seed` non-negative integer for the random generator seed.
- `parallel_evaluation` enables batch evaluations.
- `rng` user-defined Random Number Generator.

# Example

```jldoctest
julia> options = Options(f_calls_limit = 1000, debug=true, seed=1)
Options(0.0, 0.0, 0.0, 0.0, 1000.0, 0.0, 0.0, 0, false, true, true, :minimize, 0x0000000000000001)

julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> state = optimize(f, bounds, ECA(options=options))
[ Info: Initializing population...
[ Info: Starting main loop...
+=========== RESULT ==========+
| Iter.: 1
| f(x) = 6.97287
| solution.x = [-2.3628796262231875, -0.6781207370770752, -0.9642728360479853]
| f calls: 42
| Total time: 0.0004 s
+============================+

...

[ Info: Stopped since call_limit was met.
+=========== RESULT ==========+
| Iter.: 47
| f(x) = 1.56768e-08
| solution.x = [-2.2626761322304715e-5, -9.838697194048792e-5, 7.405966506272336e-5]
| f calls: 1000
| Total time: 0.0313 s
+============================+
```

"""
function Options(;
        x_tol = 0.0,
        f_tol = 0.0,
        g_tol = 0.0,
        h_tol = 0.0,
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
        seed = rand(UInt32),
        rng = Random.Xoshiro(seed)#default_rng_mh(seed),
    )


    Options(
            promote(x_tol, f_tol, g_tol, h_tol)...,
            promote(f_calls_limit, g_calls_limit, h_calls_limit, time_limit)...,
            promote(iterations)...,
            # Results options
            promote(store_convergence, show_results, debug)...,
            Symbol(search_type),
            UInt(seed),
            Bool(parallel_evaluation),
            rng
           )

end

