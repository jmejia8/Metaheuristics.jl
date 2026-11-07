include("insert_random.jl")

"""
    DE_mutation(population, F = 1.0, strategy = :rand1, best_ind = 0; rng)

Generate a `Vector` computed from population used in Differential Evolution.
Parameters: `F` is the stepsize, `strategy` can be one the following `:best1`, `:rand2`,
`:randToBest1` or `:best2`.

If the population is too small to pick distinct solutions for the strategy, some
solutions are repeated.
"""
function DE_mutation(
    population,
    F::Float64 = 1.0,
    strategy::Symbol = :rand1,
    best_ind = 0;
    rng = default_rng_mh(),
)
    x = get_position(population[1])
    v = Vector{promote_type(eltype(x), Float64)}(undef, length(x))
    DE_mutation!(v, population, F, strategy, best_ind; rng)
end

"""
    DE_mutation!(v, population, F = 1.0, strategy = :rand1, best_ind = 0; rng)

In-place version of [`DE_mutation`](@ref) that writes the mutant vector into `v`.
"""
function DE_mutation!(
    v,
    population,
    F::Float64 = 1.0,
    strategy::Symbol = :rand1,
    best_ind = 0;
    rng = default_rng_mh(),
)

    N = length(population)

    i = rand(rng, 1:N)
    # select participats
    r1 = _rand_distinct(rng, N, (i,))
    r2 = _rand_distinct(rng, N, (i, r1))

    x = get_position(population[i])
    a = get_position(population[r1])
    b = get_position(population[r2])

    # strategy is selected here
    if strategy == :rand1
        # DE/rand/1
        @. v = x + F * (a - b)
        return v
    end

    r3 = _rand_distinct(rng, N, (i, r1, r2))

    c = get_position(population[r3])

    if strategy == :rand2
        # DE/rand/2
        r4 = _rand_distinct(rng, N, (i, r1, r2, r3))
        r5 = _rand_distinct(rng, N, (i, r1, r2, r3, r4))

        d = get_position(population[r4])
        ee = get_position(population[r5])

        @. v = ee + F * (a - b + c - d)
        return v
    end

    best_ind = best_ind == 0 ? argbest(population) : best_ind
    x_best = get_position(population[best_ind])
    if strategy == :best1
        # DE/best/1
        @. v = x_best + F * (b - c)
    elseif strategy == :randToBest1
        # DE/rand-to-best/1
        @. v = x + F * (x_best - x + a - b)
    elseif strategy == :best2
        # DE/best/2
        r4 = _rand_distinct(rng, N, (i, r1, r2, r3, best_ind))
        d = get_position(population[r4])
        @. v = x_best + F * (a - b + c - d)
    else
        error("Unknown strategy " * string(strategy))
    end

    return v
end

# Random index in 1:N not in `excluded`. If every index is excluded (population too
# small), any index is returned instead of looping forever.
function _rand_distinct(rng, N, excluded::Tuple)
    _n_distinct(excluded) >= N && return rand(rng, 1:N)
    r = rand(rng, 1:N)
    while r in excluded
        r = rand(rng, 1:N)
    end
    r
end

function _n_distinct(t::Tuple)
    n = 0
    for j in eachindex(t)
        repeated = false
        for l in 1:j - 1
            t[l] == t[j] && (repeated = true; break)
        end
        repeated || (n += 1)
    end
    n
end
