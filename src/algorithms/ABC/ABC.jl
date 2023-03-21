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
  iteration: 595
    minimum: 4.03152e-28
  minimizer: [1.489845115451046e-14, 1.2207275971717747e-14, -5.671872444705246e-15]
    f calls: 30020
 total time: 0.0360 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], ABC(N = 80,  No = 20, Ne = 50, limit=5))
+=========== RESULT ==========+
  iteration: 407
    minimum: 8.94719e-08
  minimizer: [8.257485723496422e-5, 0.0002852795196258074, -3.5620824723352315e-5]
    f calls: 30039
 total time: 0.0432 s
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

    options.parallel_evaluation &&
        error("ABC is not supporting parallel evaluation. Put `options.parallel_evaluation=false`")

    D = getdim(problem)

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.iterations = parameters.N + 10000D ÷ parameters.N
        options.debug && @info "Increasing f calls limit to $(options.f_calls_limit)"
    end


    _st = gen_initial_state(problem,parameters,information,options,status)
    bees = [Bee(sol) for sol in _st.population]

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

    D = getdim(problem)
    fobj = problem.f
    bees = status.population
    Ne = parameters.Ne
    No = parameters.No
    bounds = problem.search_space
    a = bounds.lb
    b = bounds.ub

    employedPhase!(bees,problem,  Ne)
    outlookerPhase!(bees,problem, No)

    @inline genPos(D=D, a=a, b=b) = a + (b - a) .* rand(D)
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

function accuracy_stop_check(status::State{Bee}, information, options)
    cond =  !isnan(information.f_optimum) && abs(fval(status.best_sol) - information.f_optimum) < options.f_tol
    cond && (status.termination_status_code = ACCURACY_LIMIT)
    cond
end

