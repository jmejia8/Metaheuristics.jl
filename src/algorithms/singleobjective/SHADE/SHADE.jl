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
    SHADE(;N=100, information = Information(), options = Options())

Success-history based parameter adaptation for differential evolution.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> result = optimize(f, [-1 -1 -1; 1 1 1.0], SHADE())
Optimization Result
===================
  Iteration:       300
  Minimum:         3.07974e-35
  Minimizer:       [-1.73647e-18, 3.98654e-18, -3.44812e-18]
  Function calls:  30000
  Total time:      0.0647 s
  Stop reason:     Maximum objective function calls exceeded.


julia> minimizer(result)
3-element Vector{Float64}:
 -1.7364702000176986e-18
  3.9865445738182325e-18
 -3.448119784571527e-18
```


"""
function SHADE(; N::Int = 100, kargs...)
    parameters = SHADE(N, [], zeros(0), zeros(0), zeros(0), zeros(0), 1)
    Algorithm(parameters; kargs...)
end

function initialize!(
    status,
    parameters::SHADE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...,
)
    D = getdim(problem)

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug && @warn("f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    # initialize parameter history
    parameters.MCR = fill(0.5, parameters.N)
    parameters.MF = fill(0.5, parameters.N)

    return gen_initial_state(problem, parameters, information, options, status)
end

function update_state!(
    status,
    parameters::SHADE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...,
)

    new_vectors = reproduction(status, parameters, problem)

    # evaluate solutions
    new_solutions = create_solutions(new_vectors, problem, ε = options.h_tol)
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
    sort!(population, lt = (a, b) -> is_better(a, b, parameters))

    Xpb = population[[rand(max.(2, ceil.(Int, rand(N) * 0.2 * N))) for _ in 1:N]] |>
        positions
    Xr1 = population[rand(1:N, N)] |> positions
    P = eltype(population)[population; archive]
    Xr2 = P[rand(1:N, N)] |> positions

    # crossover rate
    CR = randn(N) .* sqrt(0.1) + MCR[rand(1:N, N)]
    CR = clamp.(CR, 0, 1)
    # step size
    F = randc(N) .* sqrt(0.1) + MF[rand(1:N, N)]

    for _ in 1:100
        # regenerate F when <= 0 using 100 tries (to avoid infinity loop)
        m = F .<= 0
        if any(m)
            F[m] = min.(1, randc(count(m)) .* sqrt(0.1) + MF[rand(1:N, count(m))])
        else
            break
        end
    end

    # truncate F to [0, 1] (SHADE truncates F > 1; F <= 0 left after the retries goes to 0)
    F = clamp.(F, 0, 1)

    # save the F and CR actually used, to update history in environmental_selection
    parameters.CR = CR
    parameters.F = F

    F = repeat(F, 1, D)
    mask = rand(N, D) .< CR

    # compute offspring positions (matrix X)
    X = population |> positions
    X[mask] += F[mask] .* (Xpb[mask] - X[mask] + Xr1[mask] - Xr2[mask])

    # repair solution
    for x in eachrow(X)
        reset_to_violated_bounds!(x, problem.search_space)
    end

    X
end

_shade_violation(s::AbstractConstrainedSolution) = sum_violations(s)
_shade_violation(s) = zero(fval(s))

# Improvement of the trial `child` over its `parent`, measured on the same criterion
# that made the child win under the feasibility rules: constraint violation when the
# violations differ, objective value otherwise. Positive for every successful trial.
function _shade_improvement(parent, child)
    vp, vc = _shade_violation(parent), _shade_violation(child)
    vp != vc && return vp - vc
    fval(parent) - fval(child)
end

function environmental_selection(population, parameters::SHADE)
    @assert length(population) == 2 * parameters.N

    N = parameters.N
    k = parameters.k
    MCR = parameters.MCR
    MF = parameters.MF
    CR = parameters.CR
    F = parameters.F
    new_solutions = population[N + 1:end]
    #population = population[1:N]

    survivals = Int[]
    Δ = Float64[]

    # select survival
    for (i, h) in enumerate(new_solutions)
        if is_better(h, population[i], parameters)
            push!(survivals, N + i)
            push!(Δ, _shade_improvement(population[i], h))
        else
            push!(survivals, i)
        end
    end

    if isempty(Δ)
        return survivals
    end

    mask = filter(>(N), survivals) .- N

    # update the success history only with well-defined improvements
    valid = isfinite.(Δ) .& (Δ .> 0) .& (F[mask] .> 0)
    if any(valid)
        w = Δ[valid] ./ sum(Δ[valid])
        S_CR = CR[mask[valid]]
        S_F = F[mask[valid]]
        MCR[k] = w' * S_CR
        MF[k] = (w' * S_F.^2) / (w' * S_F) # weighted Lehmer mean
        parameters.k = k % length(MCR) + 1
    end

    # update archive
    archive = parameters.archive
    append!(archive, new_solutions[mask])
    length(archive) > N && deleteat!(archive, sort(randperm(length(archive))[1:N]))

    return survivals
end
