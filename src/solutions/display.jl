import Base.show

function print_vector(io, vector)
    n_items_in_vector = 5
    if length(vector) < n_items_in_vector
        show(io, vector)
        print(io, "\n")
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
    print(io, "\n")
end

function Base.show(io::IO, solution::xf_indiv)
    @printf(io, "f = %.4e\n", solution.f)
    print(io, "x = ")

    print_vector(io, solution.x)

end


function Base.show(io::IO, solution::xfgh_indiv)
    @printf(io, "f = %.4e\n", solution.f)

    print(io, "g = ")
    print_vector(io, solution.g)
    print(io, "h = ")
    print_vector(io, solution.h)
    print(io, "x = ")
    print_vector(io, solution.x)
end

function Base.show(io::IO, solution::xFgh_indiv)
    print(io, "f = ")
    print_vector(io, solution.f)
    print(io, "g = ")
    print_vector(io, solution.g)
    print(io, "h = ")
    print_vector(io, solution.h)
    print(io, "x = ")
    print_vector(io, solution.x)

end


function Base.show(io::IO, ::MIME"text/plain", population::Array{xf_indiv})
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


function Base.show(io::IO, ::MIME"text/html", population::Array{xf_indiv})
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end


function Base.show(io::IO, ::MIME"text/plain", population::Array{xfgh_indiv})
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


function Base.show(io::IO, ::MIME"text/html", population::Array{xfgh_indiv})
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end

function Base.show(io::IO, ::MIME"text/plain", population::Array{xFgh_indiv})
    isempty(population) && show(io, population) && return

    if get(io, :compact, true)
        x = map(s -> s.f[1], population)
        y = map(s -> s.f[2], population)
        plt = scatterplot(x, y, title="F space", xlabel="f_1", ylabel="f_2")
        show(io, plt)
    else
        for sol in population
            show(io, population)
        end
    end
    return nothing
end


function Base.show(io::IO, ::MIME"text/html", population::Array{xFgh_indiv})
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end



function Base.show(io::IO, status::State)

    # if typeof(status.best_sol) != xf_indiv
    #     return println(status)
    # end

    println(io, "+=========== RESULT ==========+")
    @printf(io,"%12s %.0f\n", "iteration:", status.iteration)

    if typeof(Array(status.population)) <: Array{xFgh_indiv}
        @printf(io, "%12s", "population:")
        show(io, "text/plain", Array(status.population))

        # non-dominated
        pf = get_non_dominated_solutions(status.population)
        println(io, "\nnon-dominated solution(s):")
        show(io, "text/plain", pf)
        print(io, "\n")
    else
        @printf(io,"%12s %g\n", "minimum:", minimum(status))
        @printf(io,"%12s ", "minimizer:")
        show(io, minimizer(status))
        println(io, "")
    end



    @printf(io,"%12s %.0f\n", "f calls:", status.f_calls)
    if !isempty(status.population) &&  typeof(status.population[1]) <: Union{xfgh_indiv, xFgh_indiv}
        n = sum(map(s -> s.sum_violations ≈ 0, status.population))
        @printf(io,"%12s %d / %d in final population\n", "feasibles:", n, length(status.population))
    end
    @printf(io,"%12s %.4f s\n", "total time:", status.final_time - status.start_time)
    println(io, "+============================+")
end


function Base.show(io::IO, ::MIME"text/html", status::State)
    println(io, "<pre>")
    show(io, "text/plain", status)
    println(io, "</pre>")
end

