import Base.show

function Base.show(io::IO, solution::xf_indiv)
    @printf(io, "| f(x) = %g\n| ", solution.f)
    println(io, "| x = ", solution.x)


end

function Base.show(io::IO, solution::xfgh_indiv)
    @printf(io,"| f(x) = %g\n", solution.f)
    @printf(io,"| g(x) = ")
    println(io, solution.g)
    @printf(io,"| h(x) = ")
    println(io, solution.h)
    # @printf("| Σ max(0,g(x)) + Σ |h(x)| = ")
    # println(solution.sum_violations)
    println(io, "| x = ", solution.x)


end


function Base.show(io::IO, solution::xFgh_indiv)
    @printf(io,"| f(x) = ")
    println(io, solution.f)
    @printf(io,"| g(x) = ")
    println(io, solution.g)
    @printf(io,"| h(x) = ")
    println(io, solution.h)
    println(io, "| x = ", solution.x)


end

function Base.show(io::IO, population::Array{xFgh_indiv})
    x = map(s -> s.f[1], population)
    y = map(s -> s.f[2], population)
    plt = scatterplot(x, y, title="F space", xlabel="f_1", ylabel="f_2")
    show(io, plt)
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

