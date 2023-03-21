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



    #####################################
    # common methods
    #####################################
    # status = method.status
    information = method.information
    options = method.options
    parameters = method.parameters
    ###################################

    problem = Problem(f, search_space; parallel_evaluation=options.parallel_evaluation)
    seed!(options.seed)

    ###################################

    start_time = time()

    status = method.status
    options.debug && @info("Initializing population...")
    status = initialize!(status,parameters, problem, information, options)
    status.start_time = start_time
    method.status = status

    show_status(status, parameters, options)

    status.iteration = 1


    convergence = State{typeof(status.best_sol)}[]



    ###################################
    # store convergence
    ###################################
    if options.store_convergence
        update_convergence!(convergence, status)
    end

    options.debug && @info("Starting main loop...")

    logger(status)

    status.stop = status.stop || default_stop_check(status, information, options)
    while !status.stop
        status.iteration += 1

        update_state!(status, parameters, problem, information, options)

        # store the number of fuction evaluations
        status.f_calls = problem.f_calls

        show_status(status, parameters, options)

        if options.store_convergence
            update_convergence!(convergence, status)
        end

        status.overall_time = time() - status.start_time
        logger(status)

        # common stop criteria
        status.stop = status.stop || default_stop_check(status, information, options)

        # user defined stop criteria
        stop_criteria!(status, parameters, problem, information, options)

    end

    status.overall_time = time() - status.start_time

    final_stage!(
                 status,
                 parameters,
                 problem,
                 information,
                 options
                )

    options.debug && @info "Termination reason: " * termination_status_message(status)
    status.convergence = convergence

    return status

end


function show_status(status, parameters, options)
    !options.debug && (return)
    status.final_time = time()
    msg = "Current Status of " * string(typeof(parameters))
    @info msg
    display(status)
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

    #####################################
    # common settings
    #####################################
    information = method.information
    options = method.options
    parameters = method.parameters
    problem = Problem(f, search_space; parallel_evaluation=options.parallel_evaluation)
    ###################################

    status = method.status
    problem.f_calls = status.f_calls

    # initialization
    if status.iteration == 0 || isnothing(status.best_sol)
        options.debug && @info("Initializing population...")
        seed!(options.seed)
        start_time = time()
        status = initialize!(status,parameters, problem, information, options)
        status.start_time = start_time
        method.status = status
        status.final_time = time()
        show_status(status, parameters, options)
        return method
    end
    
    convergence = State{typeof(status.best_sol)}[]

    status.iteration += 1

    update_state!(status, parameters, problem, information, options)
    # store the number of fuction evaluations
    status.f_calls = problem.f_calls
    status.overall_time = time() - status.start_time
    status.final_time = time()

    if options.store_convergence
        update_convergence!(convergence, status)
    end

    show_status(status, parameters, options)
    logger(status)

    # common stop criteria
    status.stop = status.stop || default_stop_check(status, information, options)
    # user defined stop criteria
    stop_criteria!(status, parameters, problem, information, options)


    if status.stop
        options.debug && @info "Performing final stage due to stop criteria."
        final_stage!(
                     status,
                     parameters,
                     problem,
                     information,
                     options
                    )
    end
    

    options.debug && status.stop && @info "Should stop because: " * termination_status_message(status)

    return method
end
