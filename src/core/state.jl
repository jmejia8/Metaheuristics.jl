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
- `start_time` saves the `time()` before the optimization process.
- `final_time` saves the `time()` after the optimization process.
- `stop` if true, then stops the optimization process.

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
mutable struct State{T}
    best_sol::T
    population::Array{T}

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
    termination_status_code::TerminationStatusCode

end

function State(
        best_sol,
        population;

        # upper level parameters
        f_calls = 0,
        g_calls = 0,
        h_calls = 0,
        iteration = 1,
        success_rate = 0,
        convergence = State[],
        start_time = 0.0,
        final_time = 0.0,
        stop = false,
    )

    termination_status_code = UNKNOWN_STOP_REASON
    State(#
        best_sol,
        Array(population),

        # upper level parameters
        promote(f_calls, g_calls, h_calls, iteration)...,
        Float64(success_rate),
        State[],
        start_time,
        final_time,
        final_time - start_time,
        stop,
        termination_status_code
    )

end

"""
    minimizer(state)
Returns the approximation to the minimizer (argmin f(x)) stored in `state`.
"""
minimizer(s::State) = get_position(s.best_sol)

"""
    minimum(state::Metaheuristics.State)
Returns the approximation to the minimum (min f(x)) stored in `state`.
"""
minimum(s::State) = fval(s.best_sol)

"""
    positions(state)
If `state.population` has `N` solutions, then returns a `N`×d `Matrix`.
"""
positions(s::State) = positions(s.population)

"""
    fvals(state)
If `state.population` has `N` solutions, then returns a `Vector` with the 
objective function values from items in `state.population`.
"""
fvals(s::State{T}) where T <: AbstractSolution = fvals(s.population)

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
function convergence(s::State)
    sc = s.convergence
    isempty(sc) ? (zeros(Int, 0), zeros(0)) : (map(nfes, sc), map(minimum, sc))
end


function update_convergence!(convergence, status)
    push!(convergence, deepcopy(status))
end


termination_status_message(status::State) = termination_status_message(status.termination_status_code)
