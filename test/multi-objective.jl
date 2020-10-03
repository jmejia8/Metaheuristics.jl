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
    D = 30

    # Objective function
    f(x) = begin
        g = 1.0 + (9 / (length(x) - 1)) * sum(x[2:end])
        f1 = x[1]
        f2 = g*(1-sqrt(x[1] / g))
        return [f1, f2], [g], [0.0]
    end
    ff(x) = begin
        g = 1.0 + (9 / (length(x) - 1)) * sum(x[2:end])
        f1 = x[1]
        f2 = g*(1-sqrt(x[1] / g) - (x[1]/g) * sin(10π*x[1]))
        return [f1, f2], [g], [0.0]
    end


    bounds = Array([zeros(D) ones(D)]')

    # ECA results
    # m = ECA(N = 100, ε = 1.0, options=Options(debug=false))
    m = MOEAD_DE(D, 2, N = 300, options=Options(debug=true, iterations = 500))
    status = optimize(ff, bounds, m)
    status = optimize(ff, bounds, m)
    result = status.best_sol
    display(status)
    display(status.population)

    # print("[")
    # for s in status.population
    #     print( s.f[1], " ", s.f[2], ";")
    # end
    # println("]")

    # print("[")
    # for s in status.best_sol
    #     print( s.f[1], " ", s.f[2], ";")
    # end
    # println("]")

    @test true

end
