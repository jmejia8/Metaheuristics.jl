include("ZDT.jl")
include("DTLZ.jl")

# artificial problem
function MTP(D = 10, n_solutions = 100)
    f(x) = begin
        Q = 1.0 .+ sum((x[3:end] .- 0.8).^2)
        fx = [x[1], x[2]] .* Q
        g = [ (x[1] - 1)^2 + (x[2] - 1)^2 - 1 ]
        h = [0.0]
        return fx, g, h
    end

    bounds = Array([ zeros(D) 2.0ones(D)]')
    θ = range(π, 1.5π, length=n_solutions)

    x = 1.0 .+ cos.(θ)
    y = 1.0 .+ sin.(θ)

    pareto_set = [ generateChild(zeros(0), ([x[i], y[i]], [0.0], [0.0])) for i in 1:length(θ) ]

    return f, bounds, pareto_set


end

