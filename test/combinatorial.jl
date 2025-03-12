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
        @test isapprox(minimum(result), 0, atol=1e-3)
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

    function rkga()
        n = 5
        target_perm = collect(reverse(1:n))
        # decoder
        decode(rk) = sortperm(rk);
        # objective function
        f(rk) = sum(abs.(decode(rk) - target_perm));
        res = optimize(f, [zeros(n) ones(n)], BRKGA(num_elites=50))
        @test decode(minimizer(res)) == target_perm
    end

    function grasp_vns()
        
    end
    


    rkga()
    #### Permutation
    TSP_problem()

    #### Binary
    options = Options(seed = 1, f_tol = 1e-16, iterations=1000)
    information = Information(f_optimum = 0.0)

    ga_binary = GA(;options, information)
    simple_problem(ga_binary)

    #### Integer
    bounds = repeat([0, 10], 1, 10)
    ga_integer = GA(;mutation =PolynomialMutation(;bounds,p=1e-2),
                    crossover=SBX(;bounds),
                    environmental_selection=GenerationalReplacement(),
                    options, information
                   )

    simple_problem(ga_integer, bounds)

    #### Real
    bounds = repeat([0.0, 10], 1, 10)
    ga_real = GA(;mutation =PolynomialMutation(;bounds, p = 1e-2),
                    crossover=SBX(;bounds),
                    environmental_selection=GenerationalReplacement(),
                    options, information
                   )

    simple_problem(ga_real, bounds)
end



@testset "Combinatorial II" begin
    struct MyKPNeighborhood <: Metaheuristics.Neighborhood
        k::Int
    end

    struct KPInstance
        profit
        weight
        capacity
    end

    function Metaheuristics.compute_cost(candidates, constructor, instance::KPInstance)
        # Ration profit / weight
        ratio = instance.profit[candidates] ./ instance.weight[candidates]
        # It is assumed minimizing non-negative costs
        maximum(ratio) .- ratio
    end

    function Metaheuristics.neighborhood_structure(x, s::MyKPNeighborhood, rng)
        # this is defined due to shaking procedure requires a random one
        # not the i-th neighbor.
        i = rand(rng, 1:length(x))
        reverse!(view(x, i:min(length(x), i+s.k)))
        x
    end

    function Metaheuristics.neighborhood_structure(x, s::MyKPNeighborhood, i::Integer)
        # return the i-th neighbor around x, regarding s.k structure
        i > length(x) && return nothing
        reverse!(view(x, i:min(length(x), i+s.k)))
        x
    end

    function vns_grasp()
        profit = [55, 10,47, 5, 4, 50, 8, 61, 85, 87]
        weight = [95, 4, 60, 32, 23, 72, 80, 62, 65, 46]
        capacity = 269
        optimum = 295

        # objective function and search space
        f, search_space, _ =Metaheuristics.TestProblems.knapsack(profit, weight, capacity)
        options = Options(iterations=50, seed=1)

        ###########################################
        # VND/VNS
        ###########################################
        # list the neighborhood structures
        neighborhood_shaking = [MyKPNeighborhood(6), MyKPNeighborhood(5), MyKPNeighborhood(4)]
        neighborhood_local = [MyKPNeighborhood(3), MyKPNeighborhood(2), MyKPNeighborhood(1)]
        local_search = Metaheuristics.FirstImproveSearch()
        # instantiate VNS
        vnd = Metaheuristics.VNS(;neighborhood_shaking, neighborhood_local, local_search, options)

        res = Metaheuristics.optimize(f, search_space, vnd)
        @test -minimum(res) == optimum


        ###########################################
        # GRASP
        ###########################################
        candidates = rand(search_space)
        instance = KPInstance(profit, weight, capacity)
        constructor  = Metaheuristics.GreedyRandomizedConstructor(;candidates, instance, α = 0.95)
        local_search = Metaheuristics.BestImproveSearch()
        neighborhood = Metaheuristics.TwoOptNeighborhood()
        grasp = GRASP(;constructor, local_search)

        # optimize and display results
        res = optimize(f, search_space, grasp)
        @test -minimum(res) == optimum
    end


    vns_grasp()
end

@testset "Combinatorial: MixedInteger" begin
    function mixed_integer()
        f(sol) = sum(sol[:v]) + sum((sol[:w] .- 100).^2) - sum(sol[:x]) + sum(abs.(sol[:y] .- 1000)) + sum(sol[:z])
        # define multiple compatible search spaces
        v = BoxConstrainedSpace(-10ones(Int, 6), 10ones(Int, 6))
        w = BoxConstrainedSpace(-10ones(Int, 7), 10ones(Int, 7), rigid=false)
        x = BoxConstrainedSpace(-ones(3), ones(3))
        y = BoxConstrainedSpace(-ones(4), ones(4), rigid=false)
        z = BitArraySpace(5)

        search_space = MixedSpace(:v=>v, :w=>w, :x=>x, :y=>y, :z=>z)
        @test search_space.key_order == (:v, :w, :x, :y, :z)

        for algo in [DE(N=50), ECA(N=50), PSO(N=50), SA(), CGSA(N=50), ABC(N=50)]
            # optimize and get the results
            options = Options(verbose=false,seed=1, iterations=10)
            res = optimize(f, search_space, MixedInteger(algo; options))
            @test minimum(res) isa Number
            @test minimizer(res) isa AbstractVector
            d = Metaheuristics.vec_to_dict(minimizer(res), search_space) 
            @test d isa Dict
            @test sort(collect(keys(d))) == [:v, :w, :x, :y, :z]
        end
    end

    mixed_integer()
end
