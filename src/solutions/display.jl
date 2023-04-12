import Base.show

function print_vector(io, vector)
    n_items_in_vector = 5
    if length(vector) < n_items_in_vector
        show(io, vector)
        return
    end

    @printf(io, "[")

    for x in vector[1:2]
        @printf(io, "%1.3e, ", x)
    end
    print(io, "…, ")
    for x in vector[end:end]
        @printf(io, "%1.3e", x)
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


function describe_result(io, status::State{T}) where T <: xFgh_solution
    @printf(io, "%12s", "population:")
    show(io, "text/plain", Array(status.population))

    # non-dominated
    pf = get_non_dominated_solutions(status.population)
    println(io, "\nnon-dominated solution(s):")
    show(io, "text/plain", pf)
    print(io, "\n")
end

function describe_result(io, status::State)
    @printf(io,"%12s %g\n", "minimum:", minimum(status))
    @printf(io,"%12s ", "minimizer:")
    show(io, minimizer(status))
    println(io, "")
end


function Base.show(io::IO, status::State)
    println(io, "+=========== RESULT ==========+")
    @printf(io,"%12s %.0f\n", "iteration:", status.iteration)

    describe_result(io, status)
    @printf(io,"%12s %.0f\n", "f calls:", status.f_calls)

    if eltype(status.population) <: AbstractConstrainedSolution
        n = count(is_feasible, status.population)
        @printf(io,"%12s %d / %d in final population\n", "feasibles:", n, length(status.population))
    end
    @printf(io,"%12s %.4f s\n", "total time:", status.overall_time)

    txt = status.stop ? termination_status_message(status) : ""
    @printf(io,"%12s %s\n", "stop reason:", txt)
    println(io, "+============================+")
end


function Base.show(io::IO, ::MIME"text/html", status::State)
    println(io, "<pre>")
    show(io, "text/plain", status)
    println(io, "</pre>")
end

