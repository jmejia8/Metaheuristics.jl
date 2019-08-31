import Base.Multimedia.display

function display(solution::xf_indiv)
    @printf("| f(x) = %g\n| ", solution.f)
    @show(solution.x)


end

function display(status::State)

    if typeof(status.best_sol) != xf_indiv
        return println(status)
    end

    println("+=========== STATE ==========+")
    @printf("| Iter.: %.0f\n", status.iteration)
    display(status.best_sol)
    

    @printf("| f calls: %.0f\n", status.f_calls)
    @printf("| Total time: %.4f s\n", status.final_time - status.start_time)
    println("+============================+")
end