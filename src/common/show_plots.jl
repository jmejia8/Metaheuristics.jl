using .UnicodePlots

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
    plt = boxplot(fs, title="", xlabel="")
    show(io, plt)
    
    @printf(io, "\n%7s%14s%12s%12s", "Min.", "Mean.", "Max.", "Std.\n")
    @printf(io, "%.4e | %.4e | %.4e | %.4e\n", minimum(fs), mean(fs), maximum(fs), std(fs))

    return nothing
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


function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xf_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end




function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xfgh_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end


function Base.show(io::IO, ::MIME"text/html", population::Array{T}) where T <: xFgh_solution
    println(io, "<pre>")
    show(io, "text/plain", population)
    println(io, "</pre>")
end


Base.show(io::IO, ::MIME"text/html", pop::Vector{Bee{xf_indiv}}) = show(io, "text/html", [p.sol for p in pop])
