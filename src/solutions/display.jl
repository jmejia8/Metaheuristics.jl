import Base.show

function print_vector(io, vector)
    n_items_in_vector = 5
    if length(vector) < n_items_in_vector
        show(io, vector)
        return
    end

    @printf(io, "[")

    for x in vector[1:2]
        # show(IOContext(io, :compact=>true), x)
        # print(io, ", ")
        @printf(io, "%g, ", x)
    end
    print(io, "…, ")
    for x in vector[end:end]
        @printf(io, "%g", x)
        # show(IOContext(io, :compact=>true), x)
    end

    @printf(io, "]")
end

function Base.show(io::IO, solution::xf_indiv)
    @printf(io, "(f = %.4e, ", solution.f)
    print(io, "x = ")
    print_vector(io, solution.x)
    print(io, ")")

end


function Base.show(io::IO, solution::xfgh_indiv)
    @printf(io, "(f = %.4e", solution.f)

    print(io, ", g = ")
    print_vector(io, solution.g)
    print(io, ", h = ")
    print_vector(io, solution.h)
    print(io, ", x = ")
    print_vector(io, solution.x)
    print(io, ")")
end

function Base.show(io::IO, solution::xFgh_indiv)
    print(io, "(f = ")
    print_vector(io, solution.f)
    print(io, ", g = ")
    print_vector(io, solution.g)
    print(io, ", h = ")
    print_vector(io, solution.h)
    print(io, ", x = ")
    print_vector(io, solution.x)
    print(io, ")")

end


show_pf(io, pf) = nothing
function describe_result(io, status::State{T}) where T <: xFgh_solution
    # non-dominated
    pf = get_non_dominated_solutions(status.population)
    @printf(io,"  %-16s ", "Non-dominated:")
    println(io, length(pf), " / ", length(status.population))
    show_pf(io, pf)
end

function describe_result(io, status::State)
    @printf(io,"  %-16s %g\n", "Minimum:", minimum(status))
    @printf(io,"  %-16s ", "Minimizer:")
    print_vector(io, minimizer(status))
    println(io, "")
end


function Base.show(io::IO, status::State)
    # println(io, "+=========== RESULT ==========+")
    title = "Optimization Result"
    decor = join(repeat('=', length(title)))
    printstyled(io, title, "\n", color=:blue, bold=true)
    printstyled(io, decor, "\n", color=:gray)

    status.best_sol isa Nothing && return println(io, "  Empty status.")

    @printf(io,"  %-16s %.0f\n", "Iteration:", status.iteration)

    describe_result(io, status)
    @printf(io,"  %-16s %.0f\n", "Function calls:", status.f_calls)

    if eltype(status.population) <: AbstractConstrainedSolution
        n = count(is_feasible, status.population)
        @printf(io,"  %-16s %d / %d in final population\n", "Feasibles:", n, length(status.population))
    end
    @printf(io,"  %-16s %.4f s\n", "Total time:", status.overall_time)

    txt = status.stop ? termination_status_message(status) : ""
    @printf(io,"  %-16s %s\n", "Stop reason:", txt)
    # println(io, "+============================+")
end


function Base.show(io::IO, ::MIME"text/html", status::State)
    println(io, "<pre>")
    show(io, "text/plain", status)
    println(io, "</pre>")
end

