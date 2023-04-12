include("utils.jl")
include("before.jl")
include("during.jl")
include("after.jl")


"""
      optimize(
            f::Function, # objective function
            search_space,
            method::AbstractAlgorithm = ECA();
            logger::Function = (status) -> nothing,
      )

Minimize a n-dimensional function `f` with domain `search_space` (2×n matrix) using `method = ECA()` by default.

# Example
Minimize f(x) = Σx² where x ∈ [-10, 10]³.

Solution:

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> result = optimize(f, bounds)
+=========== RESULT ==========+
  iteration: 1429
    minimum: 2.5354499999999998e-222
  minimizer: [-1.5135301653303966e-111, 3.8688354844737692e-112, 3.082095708730726e-112]
    f calls: 29989
 total time: 0.1543 s
+============================+
```
"""
function optimize(
        f::Function, # objective function
        search_space,
        method::AbstractAlgorithm = ECA();
        logger::Function = (status) -> nothing,
    )

    status, parameters, problem, information, options, convergence = _before_optimization!(f, search_space, method, logger)

    while !status.stop
        _during_optimization!(status, parameters, problem, information, options, convergence, logger)
    end

    _after_optimization!(status, parameters, problem, information, options, convergence)

    status
end


"""
    optimize!(f, search_space, method;logger)

Perform an iteration of `method`, and save the results in `method.status`.

### Example

```julia
f, bounds, _ = Metaheuristics.TestProblems.sphere();
method = ECA()
while !Metaheuristics.should_stop(method)
    optimize!(f, bounds, method)
end
result = Metaheuristics.get_result(method)
```

See also [`optimize`](@docs).
"""
function optimize!(
        f::Function, # objective function
        search_space,
        method::AbstractAlgorithm;
        logger::Function = (status) -> nothing,
    )

    status = method.status
    # initialization
    if status.iteration == 0 || isnothing(status.best_sol)
        _before_optimization!(f, search_space, method, logger)
        return method
    end

    problem = Problem(f, search_space; parallel_evaluation=options.parallel_evaluation)
    problem.f_calls = status.f_calls

    _during_optimization!(status,
                          method.parameters,
                          problem,
                          method.information,
                          method.options,
                          empty(status.population), # convergence
                          logger)

    if status.stop
        options.debug && @info "Performing final stage due to stop criteria."
        _after_optimization!(status, parameters, problem, information, options, convergence)
    end
    

    options.debug && status.stop && @info "Should stop because: " * termination_status_message(status)

    return method
end


function optimize(
        f::Function,
        _search_space,
        ::Type{T};
        logger::Function = (status) -> nothing,
        kargs...
    ) where T <: AbstractParameters

    search_space = _mat_to_bounds(_search_space)
    # configure parameters depending on the search_space
    algo = get_parameters(f, search_space, T)
    set_user_parameters!(algo; kargs...)
    # call optimize api
    optimize(f, search_space, algo; logger)
end
