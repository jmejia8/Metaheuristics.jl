mutable struct NSGA2 <: AbstractAlgorithm
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    ε::Float64
end

"""
    NSGA2(;
        information = Information(),
        options = Options()
    )

Parameters for the metaheuristic NSGA2.

"""
function NSGA2(;
    N = 100,
    η_cr = 15,
    p_cr = 0.9,
    η_m = 20,
    p_m = -1,
    ε = eps(),
    information = Information(),
    options = Options(),
)

    parameters = NSGA2(N, promote( Float64(η_cr), p_cr, η_m, p_m, ε )...)
    Algorithm(
        parameters,
        initialize! = initialize_nsga2!,
        update_state! = update_state_nsga2!,
        is_better = is_better_nsga2,
        stop_criteria = stop_check,
        final_stage! = final_stage_nsga2!,
        information = information,
        options = options,
    )

end

function tournament_selection(P, is_better)
    a, b = rand(1:length(P)), rand(1:length(P))
    
    is_better(P[a], P[b]) || (!is_better(P[b], P[a]) && P[a].crowding > P[b].crowding ) ? P[a] : P[b]
end

function SBX_crossover(vector1, vector2, bounds, η=15, p_variable = 0.9)
    xu = view(bounds, 2,:)
    xl = view(bounds, 1,:)
    D = length(vector1)

    do_crossover = ones(Bool, D)
    do_crossover[rand(D) .> p_variable] .= false
    do_crossover[ abs.( vector1 - vector1 ) .<= eps() ] .= false

    y1 = min.( vector1, vector2 )
    y2 = max.( vector1, vector2 )
    Δ = max.(eps(), y2 - y1)

    gen_β(β) = begin
        α = 2.0 .- β .^ (-  η - 1.0 )
        R = rand(D)
        mask = R .<= 1 ./ α
        s = 1 / (η + 1)
        βq = [ mask[i] ?  (R[i] * α[i])^s : (1.0 / (2 - R[i]*α[i]))^s for i in 1:D]
        βq
    end

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β) 
    c1 = @. 0.5*(y1 + y2 -  βq*Δ)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β) 
    c2 = @. 0.5*(y1 + y2 -  βq*Δ)

    # swap
    mask = rand(Bool, D)
    cc = copy(c1)
    c1[mask] = c2[mask]
    c2[mask] = cc[mask]

    cc1 = copy(vector1)
    cc1[do_crossover] = c1[do_crossover]
    cc2 = copy(vector2)
    cc2[do_crossover] = c2[do_crossover]

    return cc1, cc2
end

# https://github.com/msu-coinlab/pymoo/blob/master/pymoo/operators/mutation/polynomial_mutation.py
function polynomial_mutation!(vector, bounds, η=20, prob = 1 / length(vector))
    do_mutation = rand(length(vector)) .< prob

    xu = view(bounds, 2,do_mutation)
    xl = view(bounds, 1,do_mutation)
    x = view(vector, do_mutation)

    δ1 = (x - xl) ./ (xu - xl)
    δ2 = (xu - x) ./ (xu - xl)

    D = length(xu)
    R = rand(D)
    mask = rand(Bool, D) .< 0.5
    s = η+1.0
    mut_pow = 1.0 / (η + 1.0)
    δq = [ mask[i] ?
            ^(2R[i] + (1 - 2R[i]) * ^(1 - δ1[i], s), mut_pow) - 1.0 :
            (2.0 * (1.0 - R[i]) + 2.0 * (R[i] - 0.5) * ^(1.0 - δ2[i], s))^mut_pow
            for i in 1:D
        ]

    vector[do_mutation] = x + δq .* ( xu - xl)
    # correct using reset to bound

end


function update_state_nsga2!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        iteration,
       )

    # fast_non_dominated_sort!(view(status.population, 1:parameters.N), engine.is_better)

    for i = 1:parameters.N

        pa = tournament_selection(status.population, engine.is_better)
        pb = tournament_selection(status.population, engine.is_better)

        # crossover
        c1, c2 = SBX_crossover( get_position(pa), get_position(pb), problem.bounds,
                              parameters.η_cr, parameters.p_cr)
       
        # mutation
        rand() < parameters.p_m && polynomial_mutation!(c1,
                                                        problem.bounds,
                                                        parameters.η_m)
        rand() < parameters.p_m && polynomial_mutation!(c2,
                                                        problem.bounds,
                                                        parameters.η_m)
       
        # rapair solutions if necesary
        replace_with_random_in_bounds!(c1, problem.bounds)
        replace_with_random_in_bounds!(c2, problem.bounds)

        # evaluate children
        child1 = generateChild(c1, problem.f(c1))
        child1.sum_violations = violationsSum(child1.g, child1.h, ε = parameters.ε)

        child2 = generateChild(c2, problem.f(c2)) 
        child2.sum_violations = violationsSum(child2.g, child2.h, ε = parameters.ε)
        status.f_calls += 2
       
        # save children
        push!(status.population, child1)
        push!(status.population, child2)
    end
    
    # non-dominated sort, crowding distance, elitist removing
    sort!(status.population, by = x -> x.rank, alg = Base.Sort.QuickSort)
    truncate_population!(status.population, parameters.N, engine.is_better)

    status.stop = engine.stop_criteria(status, information, options)
end


function initialize_nsga2!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
)
    D = size(problem.bounds, 2)

    if parameters.p_m < 0.0
        parameters.p_m = 1.0 / D
    end

    if options.iterations == 0
        options.iterations = 500
 
        if options.f_calls_limit == 0
            options.f_calls_limit = options.iterations / parameters.N + 1
        end
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        if options.iterations == 0
            options.iterations = options.f_calls_limit * parameters.N
        end
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    initialize!(problem, engine, parameters, status, information, options)

    indiv_type = typeof(status.population[1])
    if indiv_type <: xfgh_indiv || indiv_type <: xFgh_indiv
        for sol in status.population
            sol.sum_violations = violationsSum(sol.g, sol.h, ε = parameters.ε)
        end
    end


end

function final_stage_nsga2!(status, information, options)
    status.final_time = time()

    # compute Pareto front if it is a multiobjective problem
    if typeof(status.population[1].f) <: Array
        options.debug && @info "Computing Pareto front..."
        status.best_sol = get_pareto_front(status.population, is_better_nsga2)
    end
end





function is_better_nsga2(
    New::xFgh_indiv,
    Old::xFgh_indiv;
    searchType = :minimize,
    ε = 0.0,
)

    old_vio = Old.sum_violations
    new_vio = New.sum_violations

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    if searchType == :minimize
        for i in 1:length(Old.f)
            if Old.f[i] < New.f[i]
                return false
            end
        end
        return true
    end

    for i in 1:length(Old.f)
        if Old.f[i] > New.f[i]
            return false
        end
    end

    return true
end
