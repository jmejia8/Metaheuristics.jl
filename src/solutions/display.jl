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


function Base.show(io::IO, ::MIME"text/plain", population::Array{T}) where T <: xf_solution
    isempty(population) && show(io, population) && return
    
    if get(io, :compact, false)
        for sol in population
            show(io, sol)
        end
        return
    end


    print(io, "Population with ", length(population), " solutions.\n")


    fs = fvals(population)
    plt = boxplot(["f"], [fs], title="", xlabel="")
    show(io, plt)
    
    @printf(io, "\n%7s%14s%12s%12s", "Min.", "Mean.", "Max.", "Std.\n")
    @printf(io, "%.4e | %.4e | %.4e | %.4e\n", minimum(fs), mean(fs), maximum(fs), std(fs))

    return nothing
end


function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xf_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end


function Base.show(io::IO, ::MIME"text/plain", population::Array{T}) where T <: xfgh_solution
    isempty(population) && show(io, population) && return

    if get(io, :compact, false)
        for sol in population
            show(io, sol)
        end
        return
    end

    n = sum(s -> sum_violations(s) ≈ 0.0, population)
    print(io, "Population with ", n, "/", length(population)," feasible solutions.\n")

    fs = fvals(population)
    V = sum_violations.(population)
    plt = boxplot(["f", "V"], [fs, V], title="", xlabel="")
    show(io, plt)
    println(io, "\nV = Sum of contraint violations.")
end


function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xfgh_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end

function Base.show(io::IO, ::MIME"text/plain", population::Array{T}) where T <: xFgh_solution
    isempty(population) && show(io, population) && return

    if get(io, :compact, true)
        x = map(s -> s.f[1], population)
        y = map(s -> s.f[2], population)
        plt = scatterplot(x, y, title="F space", xlabel="f₁", ylabel="f₂")
        show(io, plt)
    else
        for sol in population
            show(io, population)
        end
    end
    return nothing
end


function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xFgh_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
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
    @printf(io,"%12s %.4f s\n", "total time:", status.final_time - status.start_time)

    txt = status.stop ? termination_status_message(status) : ""
    @printf(io,"%12s %s\n", "stop reason:", txt)
    println(io, "+============================+")
end


function Base.show(io::IO, ::MIME"text/html", status::State)
    println(io, "<pre>")
    show(io, "text/plain", status)
    println(io, "</pre>")
end

