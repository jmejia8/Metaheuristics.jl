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


    check_compat(search_space, method.parameters)
    # TODO: method = deepcopy(method)


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

function show_status_oneline(status, parameters, options)
    !options.verbose && (return)
    d = Any[
         "Iteration" => status.iteration,
         "Num. Evals" => nfes(status),
        ]
    # show header

    t = status.iteration
    e = nfes(status)
    m = minimum(status)
    # @printf("Iter: % 3d | Evals: % 5.0f | ", t, e)
    if m isa Number
        push!(d, "  Minimum " => m)
        # @printf("Min: % 2.3g | ", m)
    else
        n = length(get_non_dominated_solutions(status.population))
        s = sprint(print, "$n/$(length(status.population))")
        push!(d, "    NDS   " => s)
        # 
    end
    # feas = count(is_feasible.(status.population)) ÷ length(status.population)
    # @printf("Feasible: % 2d %% | ", 100feas)
    try
        push!(d, " Mean CVio" => mean(sum_violations.(status.population)))
        # @printf("CVio: % 1.2g | ", mean(sum_violations.(status.population)))
    catch
    end

    ttime = 
    push!(d, "   Time   " => @sprintf("%.4f s", status.overall_time))

    if m isa Number
        v = stop_check(status,
                       CheckConvergence(
                                        f_tol_abs = options.f_tol,
                                        f_tol_rel = options.f_tol_rel,
                                        x_tol =options.x_tol
                                       ),
                       report = false
                      ) ? "Yes" : "No"
        push!(d, "Converged " => @sprintf("%10s", v))
    end
    
    # @printf("Time:  | ", status.overall_time)

    if status.iteration <= 1 || status.iteration % 100 == 0
        nm = first.(d)
        lines = [fill('-', length(n) + 2) |> join for n in nm]
        println("+", join(lines, "+"), "+")
        println("| ", join(nm, " | "), " |")
        println("+", join(lines, "+"), "+")
    end
    print("|")

    for v in last.(d)
        if v isa Integer
            @printf("% 10d | ", v)
        elseif v isa AbstractString
            @printf("% 10s | ", v)
        elseif v isa AbstractFloat
            @printf("%1.4e | ", v)
        else 
            print(v, " | ")
        end
    end
    println("")
end

function show_status(status, parameters, options)
    !options.debug && (return show_status_oneline(status, parameters, options))
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

    # check whether optimizer is compatible with search_space
    check_compat(search_space, method.parameters)

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

"""
    get_parameters(f, search_space, H)

Get default parameters for metaheuristic `H` regarding f and the search_space.

### Example

```julia
ga = Metaheuristics.get_parameters(f, Permutations(10), GA)
```
"""
get_parameters(f, search_space, ::Type{T}) where T <: AbstractParameters = T()

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
