import Base.show

function print_vector(io, vector)
    if length(vector) < 5
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

    if get(io, :compact, true) && length(solution.x) > 5
        print_vector(io, solution.x)
    else
        show(io, solution.x)
        print(io, "\n")
    end




end


function Base.show(io::IO, solution::xfgh_indiv)
    @printf(io, "f = %.4e\n", solution.f)

    if get(io, :compact, true)
        print(io, "g = ")
        print_vector(io, solution.g)
        print(io, "h = ")
        print_vector(io, solution.h)
        print(io, "x = ")
        print_vector(io, solution.x)
    else
        show(io, solution.x)
        print(io, "\n")
    end

end

function Base.show(io::IO, solution::xFgh_indiv)

    if get(io, :compact, true)
        print(io, "f = ")
        print_vector(io, solution.f)
        print(io, "g = ")
        print_vector(io, solution.g)
        print(io, "h = ")
        print_vector(io, solution.h)
        print(io, "x = ")
        print_vector(io, solution.x)
    else
        show(io, solution.x)
        print(io, "\n")
    end

end


function Base.show(io::IO, ::MIME"text/plain", population::Array{xf_indiv})
    if get(io, :compact, false)
        for sol in population
            show(io, population)
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


function Base.show(io::IO, ::MIME"text/plain", population::Array{xfgh_indiv})

    if get(io, :compact, false)
        for sol in population
            show(io, population)
        end
        return
    end
    n = sum(s -> sum_violations(s) ≈ 0.0, population)
    print(io, length(population), "-population with ", n, " feasible solutions.\n")

    fs = fvals(population)
    V = sum_violations.(population)
    plt = boxplot(["f", "V"], [fs, V], title="", xlabel="")
    show(io, plt)
    println(io, "\nV = Sum of contraint violations.")
end

function Base.show(io::IO, ::MIME"text/plain", population::Array{xFgh_indiv})
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



function Base.show(io::IO, status::State)

    # if typeof(status.best_sol) != xf_indiv
    #     return println(status)
    # end

    println(io, "+=========== RESULT ==========+")
    @printf(io,"| Iter.: %.0f\n", status.iteration)

    if typeof(status.best_sol) <: xf_indiv || typeof(status.best_sol) <: xfgh_indiv
        show(io, status.best_sol)

    elseif typeof(Array(status.population)) <: Array{xFgh_indiv}
        println(io, "Population")
        show(io, Array(status.population))

        # non-dominated
        pf = get_non_dominated_solutions(status.population)
        println(io, "\nNon-dominated solution(s)")
        show(io, pf)
        print(io, "\n")
    end



    @printf(io,"| f calls: %.0f\n", status.f_calls)
    if !isempty(status.population) &&  typeof(status.population[1]) <: Union{xfgh_indiv, xFgh_indiv}
        n = sum(map(s -> s.sum_violations ≈ 0, status.population))
        @printf(io,"| feasibles: %d / %d in final population\n", n, length(status.population))
    end
    @printf(io,"| Total time: %.4f s\n", status.final_time - status.start_time)
    println(io, "+============================+")
end


function Base.show(io::IO, ::MIME"text/html", status::State)
    println(io, "<code>")
    show(io, "text/plain", status)
    println(io, "</code>")
end

