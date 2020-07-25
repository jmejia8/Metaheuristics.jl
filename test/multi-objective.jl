using Metaheuristics

if VERSION < v"0.7.0"
    using Base.Test
    srand(31415926534)
else
    using Test
end

# write your own tests here
@testset "Multi objective" begin


    # Dimension
    D = 5

    # Objective function
    f(x) = begin
        g = 1.0 + (9 / (length(x) - 1)) * sum(x[2:end])
        f1 = x[1]
        f2 = g*(1-sqrt(x[1] / g))
        return [f1, f2], [g], [0.0]
    end

    bounds = Array([zeros(D) ones(D)]')

    # ECA results
    status = optimize(f, bounds, ECA(ε = 1.0, options=Options(debug=false)))
    result = status.best_sol.x
    fitness = status.best_sol.f

    @test true

end
