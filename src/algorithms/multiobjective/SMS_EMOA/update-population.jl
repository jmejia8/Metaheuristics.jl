include("calc-hv.jl")

function update_population!(population, offspring, n_samples, rng = default_rng_mh())
    push!(population, offspring)

    # full sort per offspring; cheap compared with the contribution computation
    fast_non_dominated_sort!(population)

    last_front = get_last_front(population)

    # reset and update contribution in last front
    update_contribution!(population, last_front, n_samples, rng)

    deleteat!(population, argmin_contribution(population, last_front, rng))
end

# index (in population) of the solution in `last_front` with minimum contribution, ties
# broken uniformly at random (reservoir sampling)
function argmin_contribution(population, last_front, rng = default_rng_mh())
    cmin = Inf
    worst = last(last_front)
    nties = 0
    for i in last_front
        c = get_contribution(population[i])
        if c < cmin
            cmin = c
            worst = i
            nties = 1
        elseif c == cmin
            nties += 1
            rand(rng, 1:nties) == 1 && (worst = i)
        end
    end
    return worst
end

function compute_contribution(population, n_samples = 10_000, rng = default_rng_mh())
    isempty(population) && return zeros(0)

    N = length(population)
    M = length(fval(first(population)))

    # objective values, one column per solution
    P = Matrix{Float64}(undef, M, N)
    for (i, sol) in enumerate(population)
        P[:, i] .= fval(sol)
    end

    M == 2 && return hv_contrib_2d(P)
    N == 1 && return [Inf]

    # reference point: 10% beyond the nadir of the front (valid for any sign of f). All
    # solutions in a front share the same constraint violation, so feasibility is ignored.
    ref = Vector{Float64}(undef, M)
    for m in 1:M
        lo, hi = extrema(view(P, m, :))
        ref[m] = hi + (hi > lo ? 0.1 * (hi - lo) : 1.0)
    end

    M == 3 && return hv_contrib_3d(P, ref)
    return hv_contrib_mc(P, ref, n_samples, rng)
end

function update_contribution!(population, last_front, n_samples, rng = default_rng_mh())

    # reset contribution
    for sol in population
        sol.crowding = Inf
    end

    ΔS = compute_contribution(view(population, last_front), n_samples, rng)

    # save contribution for further analysis
    for (i, Δ) in zip(last_front, ΔS)
        population[i].crowding = Δ
    end

end

get_contribution(sol::xFgh_indiv) = sol.crowding

# `population` is sorted by rank: the last front is the trailing block with the same rank
function get_last_front(population)
    n = length(population)
    r = population[n].rank
    i = n
    while i > 1 && population[i - 1].rank == r
        i -= 1
    end
    return i:n
end
