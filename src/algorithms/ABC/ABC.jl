include("bee_dynamics.jl")

mutable struct ABC <: AbstractParameters
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
        information = information,
        options = options,
    )

    
end

function initialize!(
        status,
        parameters::ABC,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    
    D = size(problem.bounds, 2)

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.iterations = parameters.N + 10000D ÷ parameters.N
        options.debug && @info "Increasing f calls limit to $(options.f_calls_limit)"
    end

    bees = initialbees(parameters.N, problem)
    nevals = length(bees)

    best_sol = deepcopy(getBestBee(bees))
    population = bees
    f_calls = nevals

    return State(best_sol, population; f_calls = f_calls)
end

function update_state!(
        status,#::State{Bee},
        parameters::ABC,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    D = size(problem.bounds, 2)
    fobj = problem.f
    bees = status.population
    Ne = parameters.Ne
    No = parameters.No
    bounds = problem.bounds
    a = view(bounds, 1,:)
    b = view(bounds, 2,:)

    employedPhase!(bees,problem,  Ne)
    outlookerPhase!(bees,problem, No)

    @inline genPos(D=D, a=Array(a), b = Array(b)) = a + (b - a) .* rand(D)
    best = chooseBest(bees, status.best_sol)

    status.f_calls += Ne + No + scoutPhase!(bees, problem, genPos, parameters.limit)
    status.best_sol = best
    stop_criteria!(status, parameters, problem, information, options)

end


function final_stage!(
        status,#::State{Bee{xf_indiv}},
        parameters::ABC,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()
    # status.population = map(b -> b.sol, status.population)
    # status.best_sol = status.best_sol.sol
end

is_better_abc(bee1, bee2) = is_better(bee1.sol, bee2.sol)


function stop_criteria!(
        status,#::State{Bee},
        parameters::ABC,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    cond_budget = call_limit_stop_check(status, information, options) ||
                  iteration_stop_check(status, information, options)

    if cond_budget
        return true
    end

    cond = !isnan(information.f_optimum) &&
           abs(status.best_sol.sol.f - information.f_optimum) < options.f_tol

    options.debug && cond && @info("Stopped since accuracy was met.")
    cond
end
