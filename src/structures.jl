abstract type AbstractSolution end
abstract type AbstractAlgorithm end

mutable struct xf_indiv <: AbstractSolution # Single Objective
    x::Vector{Float64}
    f::Float64
end

mutable struct xfg_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
end

mutable struct xfgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|

end

function xfgh_indiv(
    x::Vector{Float64},
    f::Float64,
    g::Vector{Float64},
    h::Vector{Float64};
    sum_violations = 0
)
    if sum_violations <= 0
        sum_violations = violationsSum(g, h)
    end

    xfgh_indiv(x, f, g, h, sum_violations)
end

mutable struct xFgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Vector{Float64}
    g::Vector{Float64}
    h::Vector{Float64}
    rank::Int
    crowding::Float64
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
end

function xFgh_indiv(
    x::Vector{Float64},
    f::Vector{Float64},
    g::Vector{Float64},
    h::Vector{Float64};
    rank = 0,
    crowding = 0.0,
    sum_violations = 0.0
)

    if sum_violations <= 0
        sum_violations = violationsSum(g, h)
    end
    xFgh_indiv(x, f, g, h, Int(rank), crowding, sum_violations)
end

"""
    State datatype
State is used to store the current metaheuristic status. In fact, the `optimize`
function returns a `State`.

- `best_sol` Stores the best solution found so far.
- `population` is an Array{typeof(best_sol)} for population-based algorithms.
- `f_calls` is the number of objective functions evaluations.
- `g_calls`  is the number of inequality constraints evaluations.
- `h_calls` is the number of equality constraints evaluations.
- `iteration` is the current iteration.
- `success_rate` percentage of new generated solutions better that their parents. 
- `convergence` used save the `State` at each iteration.
- `start_time` saves the `time()` before the optimization proccess.
- `final_time` saves the `time()` after the optimization proccess.
- `stop` if true, then stops the optimization proccess.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> state = optimize(f, bounds)
+=========== RESULT ==========+
| Iter.: 1009
| f(x) = 7.16271e-163
| solution.x = [-7.691251412064516e-83, 1.0826961235605951e-82, -8.358428300092186e-82]
| f calls: 21190
| Total time: 0.2526 s
+============================+

julia> minimum(state)
7.162710802659093e-163

julia> minimizer(state)
3-element Array{Float64,1}:
 -7.691251412064516e-83
  1.0826961235605951e-82
 -8.358428300092186e-82
```
"""
mutable struct State
    best_sol::Any
    population::Array

    f_calls::Int
    g_calls::Int
    h_calls::Int

    iteration::Int
    success_rate::Float64
    convergence::Array{State}

    start_time::Float64
    final_time::Float64
    overall_time::Float64
    stop::Bool

end

function State(
    best_sol,
    population;

    # upper level parameters
    f_calls = 0,
    g_calls = 0,
    h_calls = 0,
    iteration = 0,
    success_rate = 0,
    convergence = State[],
    start_time = 0.0,
    final_time = 0.0,
    stop = false,
)

    State(#
        best_sol,
        Array(population),

        # upper level parameters
        promote(f_calls, g_calls, h_calls, iteration)...,
        Real(success_rate),
        State[],
        start_time,
        final_time,
        final_time - start_time,
        stop,
    )

end

"""
    minimizer(state)
Returns the approximation to the minimizer (argmin f(x)) stored in `state`.
"""
minimizer(s::State) = s.best_sol.x

"""
    minimum(state::Metaheuristics.State)
Returns the approximation to the minimum (min f(x)) stored in `state`.
"""
minimum(s::State) = s.best_sol.f

"""
    positions(state)
If `state.population` has `N` solutions, then returns a `N`×d `Matrix`.
"""
positions(s::State) = begin
    isempty(s.population) ? zeros(0,0) : Array(hcat(map(a -> a.x, s.population)...)')
end

"""
    fvals(state)
If `state.population` has `N` solutions, then returns a `Vector` with the 
objective function values from items in `state.population`.
"""
fvals(s::State) = begin
    if !isempty(s.population) && typeof(s.population[1].f) <: Vector
        return Array(hcat(map(a -> a.f, s.population)...)')
    end

    return map(a -> a.f, s.population)
end

"""
    nfes(state)
get the number of function evaluations.
"""
nfes(s::State) = s.f_calls

"""
    convergence(state)
get the data (touple with the number of function evaluations and fuction values) to plot
the convergence graph. 

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> state = optimize(f, bounds, ECA(options=Options(store_convergence=true)))
+=========== RESULT ==========+
| Iter.: 1022
| f(x) = 7.95324e-163
| solution.x = [-7.782044850211721e-82, 3.590044165897827e-82, -2.4665318114710003e-82]
| f calls: 21469
| Total time: 0.3300 s
+============================+

julia> n_fes, fxs = convergence(state);
```

"""
convergence(s::State) = begin
    sc = s.convergence
    isempty(sc) ? (zeros(Int, 0), zeros(0)) : (map(nfes, sc), map(minimum, sc))
end

mutable struct Engine
    initialize!::Function
    update_state!::Function
    is_better::Function
    stop_criteria::Function
    final_stage!::Function
end

function Engine(;
    initialize!::Function = _1(kwargs...) = nothing,
    update_state!::Function = _2(kwargs...) = nothing,
    is_better::Function = _4(kwargs...) = false,
    stop_criteria::Function = _5(kwargs...) = nothing,
    final_stage!::Function = _6(kwargs...) = nothing,
)

    Engine(initialize!, update_stat, is_better, stop_criteria, final_stage!)
end




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

    iterations::Int
    store_convergence::Bool
    show_results::Bool
    debug::Bool
    search_type::Symbol
end

"""
    Options()

`Options` stores common settings for metaheuristics such as the maximum number of iterations
debug options, maximum number of function evaluations, etc.

Main properties:

- `x_tol` tolerance to the true minimizer if specified in `Information`.
- `f_tol` tolerance to the true minimum if specified in `Information`.
- `f_calls_limit` is the maximum number of function evaluations limit.
- `iterations` is the maximum number iterationn permited.
- `store_convergence` if `true`, then push the current `State` in `State.convergence` at each generation/iteration
- `debug` if `true`, then `optimize` function reports the current `State` (and interest information) for each iterations.

# Example

```jldoctest
julia> options = Options(f_calls_limit = 1000, debug=true)
Options(0.0, 0.0, 0.0, 0.0, 1000.0, 0.0, 0.0, 0, false, true, true, :minimize)

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
| f(x) = 17.4821
| solution.x = [-0.21198067919245656, -4.0408459028606885, -1.0529741414392084]
| f calls: 42
| Total time: 0.0003 s
+============================+

...

[ Info: Stopped since call_limit was met.
+=========== RESULT ==========+
| Iter.: 47
| f(x) = 1.0924e-06
| solution.x = [-0.0008050045939313401, 0.000319255968803667, -0.0005851867103384286]
| f calls: 1000
| Total time: 0.0258 s
+============================+
```

"""
function Options(;
    x_tol::Real = 0.0,
    f_tol::Real = 0.0,
    g_tol::Real = 0.0,
    h_tol::Real = 0.0,
    f_calls_limit::Real = 0,
    g_calls_limit::Real = 0,
    h_calls_limit::Real = 0,
    iterations::Int = 0,
    store_convergence::Bool = false,
    show_results::Bool = true,
    debug::Bool = false,
    search_type::Symbol = :minimize,
)


    Options(
        promote(x_tol, f_tol, g_tol, h_tol)...,
        promote(f_calls_limit, g_calls_limit, h_calls_limit)...,
        promote(iterations)...,

        # Results options
        promote(store_convergence, show_results, debug)...,
        Symbol(search_type),
    )

end



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

julia> options = Options(f_tol = 1e-3)
Options(0.0, 0.001, 0.0, 0.0, 1000.0, 0.0, 0.0, 0, false, true, false, :minimize)

julia> state = optimize(f, bounds, ECA(information=information, options=options))
+=========== RESULT ==========+
| Iter.: 22
| f(x) = 0.000650243
| solution.x = [0.022811671589729583, 0.007052331140376011, -0.008951836265056107]
| f calls: 474
| Total time: 0.0106 s
+============================+
```
"""
function Information(;#
    f_optimum = NaN,
    x_optimum::Array{Float64} = Float64[],
)

    Information(Float64(f_optimum), x_optimum)

end


mutable struct Algorithm <: AbstractAlgorithm
    parameters::Any
    status::State
    information::Information
    options::Options
    engine::Engine
end

function Algorithm(
    parameters;
    initial_state::State = State(nothing, []),
    initialize!::Function = _1(kwargs...) = nothing,
    update_state!::Function = _2(kwargs...) = nothing,
    # is_better(a, b)  = true if x is better that y
    is_better::Function = _5(kwargs...) = false,
    stop_criteria::Function = stop_check,
    final_stage!::Function = _4(kwargs...) = nothing,
    information::Information = Information(),
    options::Options = Options(),
)


    engine = Engine(
        initialize!,
        update_state!,
        is_better,
        stop_criteria,
        final_stage!,
    )

    Algorithm(parameters, initial_state, information, options, engine)

end


struct Problem
    f::Function
    bounds::Array{Float64,2}
    g::Array{Function}
    h::Array{Function}
    type::Symbol
end

function Problem(f::Function, bounds::Array, g = Function[], h = Function[])

    type::Symbol = :constrained

    if length(g) == 0 && length(h) == 0
        type = :unconstrained
    else
        type = :constrained
    end

    Problem(f, bounds, g, h, type)
end
