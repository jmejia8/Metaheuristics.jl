"""
      optimize(
            f::Function, # objective function
            bounds::AbstractMatrix,
            method::AbstractAlgorithm = ECA();
            logger::Function = (status) -> nothing,
      )

Minimize a n-dimensional function `f` with domain `bounds` (2×n matrix) using `method = ECA()` by default.

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
        bounds::AbstractMatrix,
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

    problem = Problem(f, Array(bounds); parallel_evaluation=options.parallel_evaluation)
    seed!(options.seed)

    ###################################

    start_time = time()

    status = method.status
    options.debug && @info("Initializing population...")
    status = initialize!(status,parameters, problem, information, options)
    method.status = status

    if options.debug
        status.final_time = time()
        msg = "Current Status of " * string(typeof(parameters))
        @info msg
        display(status)
    end

    status.iteration = 1


    status.start_time = start_time
    convergence = State{typeof(status.best_sol)}[]



    ###################################
    # store convergence
    ###################################
    if options.store_convergence
        update_convergence!(convergence, status)
    end

    options.debug && @info("Starting main loop...")

    logger(status)

    while !status.stop
        status.iteration += 1

        update_state!(status, parameters, problem, information, options)

        # store the number of fuction evaluations
        status.f_calls = problem.f_calls

        if options.debug
            status.final_time = time()
            msg = "Current Status of " * string(typeof(parameters))
            @info msg
            display(status)
        end

        if options.store_convergence
            update_convergence!(convergence, status)
        end

        status.overall_time = time() - status.start_time
        logger(status)

        # common stop criteria
        status.stop = status.stop ||
        call_limit_stop_check(status, information, options) ||
        iteration_stop_check(status, information, options)  ||
        time_stop_check(status, information, options)

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

    status.convergence = convergence

    return status

end

