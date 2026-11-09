include("weights_and_ideal.jl")

mutable struct MOEAD_DE <: AbstractParameters
    nobjectives::Int
    N::Int
    F::Float64
    CR::Float64
    λ::Vector{Vector{Float64}}
    η::Float64
    p_m::Float64
    T::Int
    δ::Float64
    n_r::Float64
    z::Vector{Float64}
    B::Vector{Vector{Int}}
    s1::Float64
    s2::Float64
    τ::Float64
end

"""
    MOEAD_DE(weights ;
        F = 0.5,
        CR = 1.0,
        λ = Array{Vector{Float64}}[], # ref. points
        η = 20,
        p_m = -1.0,
        T = clamp(round(Int, 0.2*length(weights)), 2, length(weights) - 1),
        δ = 0.9,
        n_r = max(1, round(Int, 0.02*length(weights))),
        z = zeros(0),
        B = Vector{Int}[],
        s1 = 0.01,
        s2 = 20.0,
        information = Information(),
        options = Options())

`MOEAD_DE` implements the original version of MOEA/D-DE. It uses the contraint handling method
based on the sum of violations (for constrained optimizaton):
`g(x, λ, z) = max(λ .* abs.(fx - z)) + sum(max.(0, gx)) + sum(abs.(hx))`

To use MOEAD_DE, the output from the objective function should be a 3-touple
`(f::Vector, g::Vector, h::Vector)`, where `f` contains the objective functions,
`g` and `h` are the inequality and equality constraints respectively.

A feasible solution is such that `g_i(x) ≤ 0 and h_j(x) = 0`.

The number of weight vectors is the population size `N` (at least 3). `T` is the
neighborhood size (`2 ≤ T ≤ N - 1`) and `n_r ≥ 1` is the maximum number of solutions
replaced by each offspring.

Every generation creates `N` offspring. The run stops as soon as the evaluation budget
(`options.f_calls_limit`) is reached, even within a generation.

!!! note "Reusing a `MOEAD_DE` object"
    The ideal point `z` is stored in the parameters and is not reset between runs. Calling
    `optimize` again with the same object (warm start) continues from the previous
    population and the previous `z`. Create a new `MOEAD_DE` to solve a different problem.

Ref. Multiobjective Optimization Problems With Complicated Pareto Sets,
MOEA/D and NSGA-II; Hui Li and Qingfu Zhang.

# Example

Assume you want to solve the following optimizaton problem:

Minimize:

`f(x) = (x_1, x_2)`

subject to:

`g(x) = x_1^2 + x_2^2 - 1 ≤ 0`

`x_1, x_2 ∈ [-1, 1]`

A solution can be:

```julia

# Dimension
D = 2

# Objective function
f(x) = ( x, [sum(x.^2) - 1], [0.0] )

# bounds
bounds = [-1 -1;
           1  1.0
        ]

nobjectives = 2
npartitions = 100

# reference points (Das and Dennis's method)
weights = gen_ref_dirs(nobjectives, npartitions)

# define the parameters
moead_de = MOEAD_DE(weights, options=Options(debug=false, iterations = 250))

# optimize
status_moead = optimize(f, bounds, moead_de)

# show results
display(status_moead)
```

"""
function MOEAD_DE(
    weights;
    F = 0.5,
    CR = 1.0,
    λ = Array{Vector{Float64}}[],
    η = 20,
    p_m = -1.0,
    T = clamp(round(Int, 0.2 * length(weights)), 2, length(weights) - 1),
    δ = 0.9,
    n_r = max(1, round(Int, 0.02 * length(weights))),
    z::Vector{Float64} = zeros(0),
    B = Vector{Int}[],
    s1 = 0.01,
    s2 = 20.0,
    kargs...,
)

    if isempty(weights)
        error("Provide weights points")
    end

    nobjectives = length(weights[1])
    N = length(weights)

    N < 3 && throw(ArgumentError("MOEAD_DE needs at least 3 weight vectors (given $N)."))
    2 <= T <= N - 1 ||
        throw(ArgumentError("MOEAD_DE needs 2 ≤ T ≤ N - 1 = $(N - 1) (given T = $T)."))
    n_r >= 1 || throw(ArgumentError("MOEAD_DE needs n_r ≥ 1 (given n_r = $n_r)."))

    if isempty(z)
        z = fill(Inf, nobjectives)
    end

    parameters = MOEAD_DE(
        nobjectives,
        N,
        F,
        CR,
        convert(Vector{Vector{Float64}}, weights),
        η,
        p_m,
        T,
        δ,
        n_r,
        z,
        convert(Vector{Vector{Int}}, B),
        s1,
        s2,
        0.0,
    )

    initialize_closest_weight_vectors!(parameters)

    Algorithm(parameters; kargs...)

end

function initialize!(
    status,
    parameters::MOEAD_DE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...,
)

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    status = gen_initial_state(problem, parameters, information, options, status)
    n = length(status.population)
    if n != length(parameters.λ)
        throw(ArgumentError(
            "MOEAD_DE: the population size ($n) must equal the number " *
                "of weight vectors ($(length(parameters.λ))).",
        ))
    end
    parameters.nobjectives = length(status.population[1].f)

    update_reference_point!(parameters.z, status.population)
    return status

end

function update_state!(
    status::State,
    parameters::MOEAD_DE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...,
)
    # function barrier: status.population is not concretely typed
    _moead_de_generation!(
        status.population,
        status,
        parameters,
        problem,
        information,
        options,
        options.rng,
    )
    return
end

function _moead_de_generation!(
    population::AbstractVector,
    status,
    parameters::MOEAD_DE,
    problem,
    information,
    options,
    rng,
)
    N = parameters.N
    B = parameters.B
    λ = parameters.λ
    z = parameters.z

    # contiguous copies of weights, objective values and violations: the replacement
    # loop reads them for random j, which is much cheaper than chasing heap pointers
    M = length(z)
    Λ = _moead_de_columns(λ, M)
    Fpop = _moead_de_columns(map(s -> s.f, population), M)
    CV = map(sum_violations, population)

    # mating pool buffer (reused by every subproblem)
    P = Vector{Int}(undef, N)
    t0 = time()
    check_time = isfinite(options.time_limit)

    for i in 1:N
        # selection of mating pool
        if rand(rng) < parameters.δ
            Bi = B[i]
            resize!(P, length(Bi))
            copyto!(P, Bi)
        else
            resize!(P, N)
            copyto!(P, 1:N)
        end

        # reproduction
        v = MOEAD_DE_reproduction(i, P, population, parameters, problem, rng)
        # repair
        replace_with_random_in_bounds!(v, problem.search_space, rng)
        h = create_solution(v, problem)
        # function barrier: the type of h is not inferable (problem.f is a Function)
        _moead_de_update_subproblems!(population, h, P, Fpop, CV, Λ, z, parameters, rng)

        # budget-based stop within the generation; iterations are checked between them
        status.f_calls = problem.f_calls
        if call_limit_stop_check(status, information, options)
            status.stop = true
            break
        end
        if check_time && status.overall_time + (time() - t0) >= options.time_limit
            status.termination_status_code = TIME_LIMIT
            status.stop = true
            break
        end
    end
    return
end

# update z and replace up to n_r solutions of the mating pool P with the offspring h
function _moead_de_update_subproblems!(population, h, P, Fpop, CV, Λ, z, parameters, rng)
    M = length(z)
    fh = h.f
    CVh = sum_violations(h)
    update_reference_point!(z, fh)

    # threshold for the constraint handling
    τ = _moead_de_threshold(CV, P)
    parameters.τ = τ
    s1 = parameters.s1
    s2 = parameters.s2
    n_r = parameters.n_r

    c = 0
    shuffle!(rng, P)
    while c < n_r && !isempty(P)
        j = pop!(P)
        g1 = _g_column(fh, Λ, j, z)
        g2 = _g_column(Fpop, Λ, j, z)
        if _is_better_constrained_MOEAD_DE(g1, g2, CVh, CV[j], τ, s1, s2)
            population[j] = h
            @inbounds for k in 1:M
                Fpop[k, j] = fh[k]
            end
            CV[j] = CVh
            c += 1
        end
    end
    return
end

# τ = Vmin + 0.3 (Vmax - Vmin) over the mating pool, without temporaries
function _moead_de_threshold(CV::Vector{Float64}, P)
    Vmin = Vmax = CV[P[1]]
    @inbounds for k in 2:length(P)
        V = CV[P[k]]
        Vmin = min(Vmin, V)
        Vmax = max(Vmax, V)
    end
    return Vmin + 0.3 * (Vmax - Vmin)
end

# M×N matrix whose columns are the given vectors
function _moead_de_columns(vs, M)
    A = Matrix{Float64}(undef, M, length(vs))
    for (j, v) in enumerate(vs)
        length(v) == M || throw(DimensionMismatch("expected vectors of length $M"))
        @inbounds for k in 1:M
            A[k, j] = v[k]
        end
    end
    A
end

"""
    MOEAD_DE_reproduction(a, b, c, F, CR, p_m, η, bounds[, rng])

Perform Differential Evolution operators and polynomial mutation using three vectors
`a, b, c` and parameters `F, CR, p_m, η`, i.e., stepsize, crossover and
mutation probability.
"""
function MOEAD_DE_reproduction(
    a,
    b,
    c,
    F,
    CR,
    p_m,
    η,
    bounds::BoxConstrainedSpace,
    rng = default_rng_mh(),
)
    D = length(a)
    v = zeros(D)

    lb = bounds.lb
    ub = bounds.ub
    e = 1.0 / (η + 1)

    @inbounds for j in 1:D
        # binomial crossover
        if rand(rng) < CR
            v[j] = a[j] + F * (b[j] - c[j])
        else
            v[j] = a[j]
        end
        # polynomial mutation
        if rand(rng) < p_m
            r = rand(rng)
            if r < 0.5
                σ_k = (2.0 * r)^e - 1
            else
                σ_k = 1 - (2.0 - 2.0 * r)^e
            end
            v[j] = v[j] + σ_k * (ub[j] - lb[j])
        end
    end

    v
end

function MOEAD_DE_reproduction(
    i,
    P_idx,
    population,
    parameters::MOEAD_DE,
    problem,
    rng = default_rng_mh(),
)
    # select participants
    r1 = i
    r2 = rand(rng, P_idx)
    while r1 == r2
        r2 = rand(rng, P_idx)
    end

    r3 = rand(rng, P_idx)
    while r3 == r1 || r3 == r2
        r3 = rand(rng, P_idx)
    end

    a = get_position(population[r1])
    b = get_position(population[r2])
    c = get_position(population[r3])

    MOEAD_DE_reproduction(
        a,
        b,
        c,
        parameters.F,
        parameters.CR,
        parameters.p_m < 0.0 ? 1.0 / getdim(problem.search_space) : parameters.p_m,
        parameters.η,
        problem.search_space,
        rng,
    )
end

g_te_ap(gx, V, τ, s1, s2) = V < τ ? gx + s1 * V^2 : gx + s1 * τ^2 + s2 * (V - τ)

function is_better_constrained_MOEAD_DE(g1, g2, sol1, sol2, parameters)
    _is_better_constrained_MOEAD_DE(
        g1,
        g2,
        sum_violations(sol1),
        sum_violations(sol2),
        parameters.τ,
        parameters.s1,
        parameters.s2,
    )
end

function _is_better_constrained_MOEAD_DE(g1, g2, CV1, CV2, τ, s1, s2)
    if CV1 == 0.0 && CV2 == 0.0
        return g1 <= g2
    end
    return g_te_ap(g1, CV1, τ, s1, s2) <= g_te_ap(g2, CV2, τ, s1, s2)
end

function final_stage!(
    status::State,
    parameters::MOEAD_DE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...,
)

    status.final_time = time()
end
