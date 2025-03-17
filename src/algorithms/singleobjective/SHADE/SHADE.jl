mutable struct SHADE <: AbstractDifferentialEvolution
    N::Int
    archive::Vector
    MF::Vector{Float64}
    F::Vector{Float64}
    MCR::Vector{Float64}
    CR::Vector{Float64}
    k::Int
end


"""

"""
function SHADE(; N::Int = 100, kargs...)
    parameters = SHADE(N, [], zeros(0), zeros(0),zeros(0), zeros(0), 1)
    Algorithm(parameters; kargs...)
end

function initialize!(
        status,
        parameters::SHADE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    D = getdim(problem)


    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
        @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    # initialize parameter history
    parameters.MCR = fill(0.5, parameters.N)
    parameters.MF = fill(0.5, parameters.N)

    return gen_initial_state(problem,parameters,information,options,status)
end


function update_state!(
        status,
        parameters::SHADE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    new_vectors = reproduction(status, parameters, problem)

    # evaluate solutions
    new_solutions = create_solutions(new_vectors, problem,ε=options.h_tol)
    append!(status.population, new_solutions)

    environmental_selection!(status.population, parameters)
    status.best_sol = get_best(status.population)
end


function reproduction(status, parameters::SHADE, problem)
    population = status.population
    @assert !isempty(population)
    archive = parameters.archive
    MCR = parameters.MCR
    MF = parameters.MF

    # Cauchy distribution generator
    randc(n) = randn(n) ./ randn(n) # X/Y ~ Cauchy(0, 1)

    # population size
    N = parameters.N
    # number of decisions
    D = getdim(problem)
    sort!(population, lt=(a, b) -> is_better(a, b, parameters))

    Xpb = population[[rand(max.(2, ceil.(Int, rand(N)*0.2*N))) for _ = 1:N]] |> positions
    Xr1 = population[rand(1:N, N)] |> positions
    P   = eltype(population)[population ; archive]
    Xr2 = P[rand(1:N, N)] |> positions

    # crossover rate
    CR  = randn(N).*sqrt(0.1) + MCR[rand(1:N, N)]
    CR  = clamp.(CR, 0, 1)
    # step size
    F = randc(N).*sqrt(0.1) + MF[rand(1:N, N)]

    for _ in 1:100 # regenerate F when <= 0 using 100 tries (to avoid infinity loop)
        m = F .<= 0
        if any(m)
            F[m] = min.(1, randc(count(m)).*sqrt(0.1) + MF[rand(1:N, count(m))])
        else
            break
        end
    end

    # save F and CR to update history in environmental_selection
    parameters.CR = CR
    parameters.F = F

    F = repeat(clamp.(F, 0, 1), 1, D)
    mask = rand(N, D) .< CR

    # compute offspring positions (matrix X)
    X = population |> positions
    X[mask] += F[mask] .* (Xpb[mask]-X[mask]+Xr1[mask]-Xr2[mask])

    X
end


function environmental_selection(population, parameters::SHADE)
    @assert  length(population) == 2*parameters.N

    N = parameters.N
    k = parameters.k
    MCR = parameters.MCR
    MF = parameters.MF
    CR = parameters.CR
    F = parameters.F
    new_solutions = population[N+1:end]
    #population = population[1:N]

    survivals = Int[]
    Δ = Float64[]

    # select survival
    for (i, h) in enumerate(new_solutions)
        if is_better(h, population[i], parameters)
            push!(survivals, N + i)
            push!(Δ, fval(population[i]) - fval(h))
        else
            push!(survivals, i)
        end
    end

    if isempty(Δ)
        return survivals
    end

    mask = filter(>(N), survivals) .- N

    w      = Δ ./ sum(Δ)
    MCR[k] = w'*CR[mask]
    MF[k]  = (w'*F[mask].^2) / (w'*F[mask])
    k      = k % length(MCR) + 1

    # update archive
    archive  = parameters.archive
    append!(archive, new_solutions[mask])
    length(archive) > N && deleteat!(archive, sort(randperm(length(archive))[1:N]))
    

    return survivals
end

