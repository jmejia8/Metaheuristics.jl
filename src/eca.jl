
mutable struct ECA <: AbstractAlgorithm
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

function ECA(;
    η_max::Float64 = 2.0,
    K::Int = 7,
    N::Int = 0,
    N_init::Int = N,
    p_exploit::Float64 = 0.95,
    p_bin::Float64 = 0.02,
    ε::Float64 = 0.0,
    p_cr::Array{Float64} = Float64[],
    adaptive::Bool = false,
    resize_population::Bool = false,
    information = Information(),
    options = Options(),
)


    N_init = N


    parameters = ECA(
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
        initialize! = initialize_eca!,
        update_state! = update_state_eca!,
        is_better = is_better_eca,
        stop_criteria = stop_check,
        final_stage! = final_stage_eca!,
        information = information,
        options = options,
    )

end

function fitnessToMass(fitness::Vector{Float64}, searchType::Symbol)
    m = minimum(fitness)

    if m < 0
        fitness = 2 * abs(m) .+ fitness
    end

    if searchType == :minimize
        fitness = 2 * maximum(fitness) .- fitness
    end

    return fitness
end

function getMass(U::Array{xf_indiv,1}, searchType::Symbol; kargs...)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness, searchType)
end

function getMass(U::Array{xfg_indiv,1}, searchType; kargs...)
    return getMass(U, searchType)

end

function getMass(U::Array{xfgh_indiv,1}, searchType; ε = 0.0)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        v = U[i].sum_violations
        if v > 0.0
            fitness[i] = v
        else
            fitness[i] = U[i].f
        end
    end

    return fitnessToMass(fitness, searchType)

end

function getMass(U::Array{xFgh_indiv,1}, searchType; ε = 0.0)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        v = U[i].sum_violations
        if v > 0.0
            fitness[i] = v
        else
            fitness[i] = sum(U[i].f)
        end
    end

    return fitnessToMass(fitness, searchType)

end

function center(U, mass)
    d = length(U[1].x)

    c = zeros(Float64, d)

    for i = 1:length(mass)
        c += mass[i] .* U[i].x
    end

    return c / sum(mass)
end

function center(U::Array, searchType::Symbol; ε = 0.0)
    n = length(U)

    mass = getMass(U, searchType; ε = ε)

    return center(U, mass),
    getWorstInd(U, searchType, is_better_eca),
    getBestInd(U, searchType, is_better_eca)
end

function getU(P::Array, K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N - K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j.+1]
    end

    return P[U_ids]
end

function crossover(
    x::Vector{Float64},
    y::Vector{Float64},
    p_cr::Vector{Float64},
)
    D = length(x)
    tmp2 = zeros(D)
    for j = 1:D
        if rand() < p_cr[j]
            y[j] = x[j]
            tmp2[j] += 1
        end
    end

    return y, tmp2
end

function adaptCrossover(p_cr::Vector{Float64}, M::Vector{Float64})
    p_best = p_cr[indmin(M)]

    for i = 1:length(p_cr)
        if M[i] > 0.3
            pn = abs(p_best .+ 0.3 * randn())
            if pn > 1.0
                pn = 1.0
            end

            if pn < 0.0
                pn = p_best
            end
            p_cr[i] = pn
        end
    end

    return p_cr
end

function resizePop!(P::Array, N_new::Int, K::Int; is_better = is_better)
    N = length(P)

    if N == N_new
        return P
    end

    f = zeros(N)
    for i = 1:N
        f[i] = P[i].f
    end

    ids = sortperm(P, lt = is_better)[1:N_new]
    return P[ids]
end

function update_state_eca!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    iteration,
)
    K = parameters.K
    I = randperm(parameters.N)
    D = size(problem.bounds, 2)
    is_multiobjective = typeof(status.population[1]) == xFgh_indiv

    parameters.adaptive && (Mcr_fail = zeros(D))

    stop = false
    a = problem.bounds[1, :]
    b = problem.bounds[2, :]
    ε = 0.0

    if is_multiobjective
        shuffle!(status.population)
    end

    # For each elements in Population
    for i = 1:parameters.N
        p = status.f_calls / options.f_calls_limit

        if parameters.ε > 0
            ε = parameters.ε * (p - 1)^4

        else
            ε = 0.0
        end

        # current
        x = status.population[i].x

        # generate U masses
        U = getU(status.population, parameters.K, I, i, parameters.N)

        # generate center of mass
        c, u_worst, u_best = center(U, :minimize)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        if p < parameters.p_exploit
            # u: worst element in U
            u = U[u_worst].x

            # current-to-center/bin
            y = x .+ η .* (c .- u)
        elseif parameters.p_exploit < 0
            y =
                x .+ (1 - p^5) * η * (c .- u) .+
                (p^5) * η * (status.best_sol.x .- c)
        else
            # current-to-best/bin
            y = x .+ η .* (status.best_sol.x .- c)
        end

        # binary crossover
        y, M_current = crossover(U[u_best].x, y, parameters.p_cr)

        y = correct(y, c, a, b)

        sol = generateChild(y, problem.f(y))
        if ε > 0.0
            sol.sum_violations = violationsSum(sol.g, sol.h, ε = ε)
        end

        status.f_calls += 1

        # replace worst element
        if engine.is_better(sol, status.population[i])
            wi = getWorstInd(status.population, :minimize, (w,z) -> engine.is_better(w, z))
            status.population[wi] = sol
            if engine.is_better(sol, status.best_sol)
                status.best_sol = sol
            end
        else
            if is_multiobjective && !engine.is_better(status.population[i], sol)
                push!(status.population, sol)
            end
            parameters.adaptive && (Mcr_fail += M_current)
        end

        status.stop = engine.stop_criteria(status, information, options)
        status.stop && break
    end

    if status.stop
        return
    end

    if parameters.adaptive
        parameters.p_cr =
            adaptCrossover(parameters.p_cr, Mcr_fail / parameters.N)
    end

    if is_multiobjective && length(status.population) > parameters.N
        status.population = truncate_population!(status.population, parameters.N, (w,z) -> engine.is_better(w, z))
    end

    p = status.f_calls / options.f_calls_limit

    if parameters.resize_population
        K = parameters.K

        # new size
        parameters.N = 2K .+ round(Int, (1 - p) * (parameters.N_init .- 2K))

        if parameters.N < 2K
            parameters.N = 2K
        end

        status.population = resizePop!(
            status.population,
            parameters.N,
            K,
            (a, b) -> engine.is_better(a, b),
        )
    end
end


function initialize_eca!(
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

    N_init = parameters.N


    if parameters.adaptive
        parameters.p_cr = rand(D)
    else
        parameters.p_cr = parameters.p_bin .* ones(D)
    end

end

function final_stage_eca!(status, information, options)
    status.final_time = time()

    # compute Pareto front if it is a multiobjective problem
    if typeof(status.population[1].f) <: Array
        options.debug && @info "Computing Pareto front..."
        status.best_sol = get_pareto_front(status.population, is_better_eca)
    end
end


is_better_eca(
    New::xf_indiv,
    Old::xf_indiv;
    searchType = :minimize,
    leq = false,
    kargs...,
) = New.f < Old.f


function is_better_eca(
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


function is_better_eca(
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
