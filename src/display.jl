import Base.Multimedia.display

function display(solution::xf_indiv)
    @printf("| f(x) = %g\n| ", solution.f)
    @show(solution.x)


end

function display(solution::xfgh_indiv)
    @printf("| f(x) = %g\n", solution.f)
    @printf("| g(x) = ")
    println(solution.g)
    @printf("| h(x) = ")
    println(solution.h)
    @show(solution.x)


end

function display(population::Array{xFgh_indiv})
    x = map(s -> s.f[1], population)
    y = map(s -> s.f[2], population)
    plt = scatterplot(x, y, title="Population", xlabel="f_1", ylabel="f_2")
    display(plt)
end

function display(status::State)

    # if typeof(status.best_sol) != xf_indiv
    #     return println(status)
    # end

    println("+=========== STATE ==========+")
    @printf("| Iter.: %.0f\n", status.iteration)
    display(status.population)


    @printf("\n| f calls: %.0f\n", status.f_calls)
    @printf("| Total time: %.4f s\n", status.final_time - status.start_time)
    println("+============================+")
end
