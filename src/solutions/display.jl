import Base.Multimedia.display

function display(solution::xf_indiv)
    @printf("| f(x) = %g\n| ", solution.f)
    println("| x = ", solution.x)


end

function display(solution::xfgh_indiv)
    @printf("| f(x) = %g\n", solution.f)
    @printf("| g(x) = ")
    println(solution.g)
    @printf("| h(x) = ")
    println(solution.h)
    println("| x = ", solution.x)


end


function display(solution::xFgh_indiv)
    @printf("| f(x) = ")
    println(solution.f)
    @printf("| g(x) = ")
    println(solution.g)
    @printf("| h(x) = ")
    println(solution.h)
    println("| x = ", solution.x)


end

function display(population::Array{xFgh_indiv})
    x = map(s -> s.f[1], population)
    y = map(s -> s.f[2], population)
    plt = scatterplot(x, y, title="F space", xlabel="f_1", ylabel="f_2")
    display(plt)
end

function display(status::State)

    # if typeof(status.best_sol) != xf_indiv
    #     return println(status)
    # end

    println("+=========== RESULT ==========+")
    @printf("| Iter.: %.0f\n", status.iteration)

    if typeof(status.best_sol) <: xf_indiv || typeof(status.best_sol) <: xfgh_indiv
        display(status.best_sol)

    elseif typeof(Array(status.population)) <: Array{xFgh_indiv}
        println("Population")
        display(Array(status.population))

        b = nothing
        try
            b = Array(status.best_sol)
        catch
            b = nothing
        end
        if !isnothing(b) && !isempty(b)
            println("\nBest solution(s)")
            display(b)
        end
        print("\n")
    end



    @printf("| f calls: %.0f\n", status.f_calls)
    if typeof(status.population[1]) <: Union{xfgh_indiv, xFgh_indiv}
        n = sum(map(s -> s.sum_violations ≈ 0, status.population))
        @printf("| feasibles: %d / %d in final population\n", n, length(status.population))
    end
    @printf("| Total time: %.4f s\n", status.final_time - status.start_time)
    println("+============================+")
end
