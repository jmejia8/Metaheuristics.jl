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

julia> state = optimize(f, bounds, ECA(options=Options(seed=1)))
Optimization Result
===================
  Iteration:       303
  Minimum:         7.02515e-49
  Minimizer:       [7.93381e-25, -2.38145e-25, 1.27862e-25]
  Function calls:  6363
  Total time:      0.0349 s
  Stop reason:     Due to Convergence Termination criterion.

julia> minimum(state)
7.025152289113865e-49

julia> minimizer(state)
3-element Vector{Float64}:
  7.933813096004835e-25
 -2.381454635702266e-25
  1.278624129134673e-25
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

julia> state = optimize(f, bounds, ECA(options=Options(seed=1, store_convergence=true)))
Optimization Result
===================
  Iteration:       303
  Minimum:         7.02515e-49
  Minimizer:       [7.93381e-25, -2.38145e-25, 1.27862e-25]
  Function calls:  6363
  Total time:      0.0349 s
  Stop reason:     Due to Convergence Termination criterion.

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
