"""
    get_parameters(f, search_space, H)

Get default parameters for metaheuristic `H` regarding f and the search_space.

### Example

```julia
ga = Metaheuristics.get_parameters(f, Permutations(10), GA)
```
"""
get_parameters(f, search_space, ::Type{T}) where T <: AbstractParameters = T()


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
    if m isa Number
        push!(d, "  Minimum " => m)
    else
        n = length(get_non_dominated_solutions(status.population))
        s = sprint(print, "$n/$(length(status.population))")
        push!(d, "    NDS   " => s)
        # 
    end
    # feas = count(is_feasible.(status.population)) ÷ length(status.population)
    try
        push!(d, " Mean CVio" => mean(sum_violations.(status.population)))
    catch
    end

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
    

    if status.iteration <= 1 || status.iteration % 1000 == 0
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
