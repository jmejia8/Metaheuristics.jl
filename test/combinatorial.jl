@testset "Combinatorial" begin

    function test_results(result)
        show(IOBuffer(), "text/html", result)
        show(IOBuffer(), "text/plain", result.population)
        show(IOBuffer(), "text/html", result.population)
        show(IOBuffer(), result.population[1])
    end

    function simple_problem(method, bounds = [zeros(Bool,10) ones(Bool, 10)])
        ff(x) = sum( abs.(x .- 1.0) )

        # number of function evaluations
        f_calls = 0

        f(x) = begin
            f_calls += 1
            ff(x)
        end

        result = optimize(f, bounds, method)
        @test minimum(result) ≈ 0
        test_results(result)
        result
    end

    function TSP_problem()
        # TST Problem
        θ = range(0, 2π, length = 10)
        n = length(θ)
        coor = [cos.(θ) sin.(θ)]
        _dist(x, y) = sum( (x-y).^2 )
        Dists = [_dist(coor[i,:], coor[j,:]) for i in 1:n, j in 1:n]
        f(path) = sum( Dists[path[i],path[i+1]] for i in 1:n-1 ) + Dists[path[end],path[1]]

        bounds = [fill(1, n) fill(n, n)]

        # general options
        options = Options(iterations = 100, f_calls_limit=1e10, seed=1, f_tol=1e-8)
        information = Information(f_optimum = f(1:n))

        # parameters
        ga = GA(;
                initializer = RandomPermutation(N = 100),
                crossover = OrderCrossover(),
                mutation = SlightMutation(),
                options,
                information
               )

        result = optimize(f, bounds, ga)
        @test minimum(result) ≈ f(1:n)
        test_results(result)
    end
 
    #### Permutation
    TSP_problem()
 
    #### Binary
    options = Options(seed = 1, f_tol = 1e-16)
    information = Information(f_optimum = 0.0)

    ga_binary = GA(;options, information)
    simple_problem(ga_binary)

    #### Integer
    bounds = repeat([0, 10], 1, 10)
    ga_integer = GA(;mutation =PolynomialMutation(;bounds),
                    crossover=SBX(;bounds),
                    environmental_selection=GenerationalReplacement(),
                    options, information
                   )

    simple_problem(ga_integer, bounds)
end


