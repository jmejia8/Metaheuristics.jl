include("bee_dynamics.jl")

mutable struct ABC
    N::Int
    Ne::Int
    No::Int
    limit::Int
end

"""
    ABC(;
        N = 50,
        Ne = div(N+1, 2),
        No = div(N+1, 2),
        limit=10,
        information = Information(),
        options = Options()
    )

ABC implements the original parameters for the Artificial Bee Colony Algorithm.
`N` is the population size, `Ne` is the number of employees, `No` is the
number of outlookers bees. `limit` is related to the times that a solution
is visited.


# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], ABC())
+=========== RESULT ==========+
| Iter.: 593
| f(x) = 3.54833e-25
| solution.x = [3.448700205761237e-13, 4.805851037329074e-13, 7.025504722610375e-14]
| f calls: 30019
| Total time: 0.2323 s
+============================+
julia> optimize(f, [-1 -1 -1; 1 1 1.0], ABC(N = 80,  No = 20, Ne = 50, limit=5))
+=========== RESULT ==========+
| Iter.: 405
| f(x) = 2.24846e-07
| solution.x = [0.0002682351072804559, 0.00020460896416511776, 0.0003332131896109299]
| f calls: 30043
| Total time: 0.2652 s
+============================+
```

"""
function ABC(;
        N = 50,
        Ne = div(N+1, 2),
        No = div(N+1, 2),
        limit=10,
        information = Information(),
        options = Options()
    )
    
    parameters = ABC(N, Ne, No, limit)

    Algorithm(
        parameters,
        initialize! = initialize_abc!,
        update_state! = update_state_abc!,
        is_better = is_better_abc,
        stop_criteria = stop_check_abc,
        final_stage! = final_stage_abc!,
        information = information,
        options = options,
    )

    
end

function initialize_abc!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
   )
    
    D = size(problem.bounds, 2)

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.iterations = parameters.N + 10000D ÷ parameters.N
        options.debug && @info "Increasing f calls limit to $(options.f_calls_limit)"
    end

    bees = initialbees(problem.f, parameters.N, problem.bounds)
    nevals = length(bees)

    status.best_sol = deepcopy(getBest(bees))
    status.population = bees
    status.f_calls = nevals
end

function update_state_abc!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        iteration,
       )

    D = size(problem.bounds, 2)
    fobj = problem.f
    bees = status.population
    Ne = parameters.Ne
    No = parameters.No
    bounds = problem.bounds
    a = view(bounds, 1,:)
    b = view(bounds, 2,:)

    employedPhase!(bees, fobj, Ne, bounds)        
    outlookerPhase!(bees, fobj, No, bounds)

    @inline genPos(D=D, a=Array(a), b = Array(b)) = initializeSol(D, a, b)
    best = chooseBest(bees, status.best_sol)

    status.f_calls += Ne + No + scoutPhase!(bees, fobj, genPos, parameters.limit)
    status.best_sol = best
    status.stop = engine.stop_criteria(status, information, options)

end


function final_stage_abc!(status, information, options)
    status.final_time = time()
    status.population = map(b -> b.sol, status.population)
    # status.best_sol = status.best_sol.sol
end

is_better_abc(bee1, bee2) = is_better(bee1.sol, bee2.sol)


function stop_check_abc(status, information, options)
    cond2 = call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)

    cond = cond2 || !isnan(information.f_optimum) && abs(status.best_sol.f - information.f_optimum) < options.f_tol

    options.debug && cond && @info("Stopped since accuracy was met.")
    cond
end
