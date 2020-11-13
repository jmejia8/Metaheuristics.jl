mutable struct NSGA2 <: AbstractAlgorithm
    η_max::Float64
    K::Int
    N::Int
    N_init::Int
    p_exploit::Float64
    p_bin::Float64
    ε::Float64
    p_cr::Array{Float64}
    adaptive::Bool
    resize_population::Bool
end

"""
    NSGA2(;
        information = Information(),
        options = Options()
    )

Parameters for the metaheuristic NSGA2.

"""
function NSGA2(;

    information = Information(),
    options = Options(),
)


    N_init = N


    parameters = NSGA2(
        η_max,
        K,
        N,
        N_init,
        p_exploit,
        p_bin,
        ε,
        p_cr,
        adaptive,
        resize_population,
    )
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

function tournament_selection(P)
    a, b = rand(1:length(P)÷2), rand(1:length(P)÷2)
    P[a] < P[b] ? P[a] : P[b]
end

function SBX_crossover(vector1, vector2, η, prob_per_variable = 0.5)
    xu = view(bounds, 2,:)
    xl = view(bounds, 1,:)
    D = length(vector1)

    do_crossover = rand(Bool, D)
    do_crossover[ abs.( vector1 - vector1 ) .<= eps() ] .= false

    y1 = min.( vector1, vector2 )
    y2 = max.( vector1, vector2 )
    Δ = max.(eps(), y2 - y1)

    gen_β(β) = begin
        α = 2.0 - β .^ (-  η - 1.0 )
        R = rand(D)
        mask = R <= 1 ./ α
        s = 1 / (η + 1)
        βq = [ mask[i] ?  (R[i] * α[i])^s : (1.0 / (2 - R[i]*α[i]))^s for i in 1:D]
        βq
    end

    β = 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β) 
    c1 = @. 0.5*(y + y2 -  βq*Δ)

    β = 1.0 + (2.0 * (y1 - xl) / δ)
    βq = gen_β(β) 
    c2 = @. 0.5*(y + y2 -  βq*Δ)

    # swap
    mask = rand(Bool, D)
    cc = copy(c1)
    c1[mask] = c2[mask]
    c2[mask] = cc
    return c1, c2
end

# https://github.com/msu-coinlab/pymoo/blob/master/pymoo/operators/mutation/polynomial_mutation.py
function polynomial_mutation!(vector, bounds)
    xu = view(bounds, 2,:)
    xl = view(bounds, 1,:)

    δ1 = (vector - xl) / (xu - xl)
    δ2 = (xu - vector) / (xu - xl)

    idx = rand(length(vector)) <= 0.5

    fill!(vector, vector + δ .* ( xu - xl))

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

    fast_non_dominated_sort!(view(status.population, 1:parameters.N))

    for i = 1:2:parameters.N

        pa = tournament_selection(P)
        pb = tournament_selection(P)

        c1, c2 = SBX_crossover( get_position(pa), get_position(pb) )
        
        rand() < parameters.p_m && mutate!(c1)
        rand() < parameters.p_m && mutate!(c2)

        eval!(P[popSize + i], z, fCV)
        eval!(P[popSize + i + 1], z, fCV)
    end

    fast_non_dominated_sort!(status.population)
    sort!(status.population, by = x -> x.rank, alg = Base.Sort.QuickSort)


    status.population = truncate_population!(status.population, parameters.N, (w,z) -> engine.is_better(w, z))
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


    if parameters.N <= parameters.K
        parameters.N = parameters.K * D
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
            @warn( "f_calls_limit increased to $(options.f_calls_limit)")
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
    N_init = parameters.N


    if parameters.adaptive
        parameters.p_cr = rand(D)
    else
        parameters.p_cr = parameters.p_bin .* ones(D)
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


is_better_nsga2(
    New::xf_indiv,
    Old::xf_indiv;
    searchType = :minimize,
    leq = false,
    kargs...,
) = New.f < Old.f


function is_better_nsga2(
    New::xfgh_indiv,
    Old::xfgh_indiv;
    searchType = :minimize
)

    old_vio = Old.sum_violations
    new_vio = New.sum_violations

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    if searchType == :minimize
        return New.f < Old.f
    end


    return New.f > Old.f
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
