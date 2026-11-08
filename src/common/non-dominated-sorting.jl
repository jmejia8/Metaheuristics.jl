"""
    non_dominated_sort(population)

Return a vector of integers `r` containing in `r[i]` the rank for `population[i]`.
"""
function non_dominated_sort(pop)
    _ens_applicable(pop) && return _non_dominated_sort_ens(pop)
    return _non_dominated_sort_naive(pop)
end

# ENS-BS (Zhang et al. 2015, doi:10.1109/TEVC.2014.2308305). Solutions are visited in
# an order compatible with (constrained) dominance, hence a solution can only be
# dominated by previously visited ones. Requires that dominance is a strict partial
# order (not guaranteed when NaN values are present).
_ens_applicable(pop) = false
_ens_applicable(pop::AbstractVector{<:Vector{<:Real}}) = all(v -> !any(isnan, v), pop)
function _ens_applicable(pop::AbstractVector{<:xFgh_solution})
    all(s -> !isnan(s.sum_violations) && !any(isnan, s.f), pop)
end

# if `a` (constrained-)dominates `b` then `_lex_less(a, b)`
@inline function _lex_less(a::AbstractVector, b::AbstractVector)
    @inbounds for i in eachindex(a, b)
        a[i] < b[i] && return true
        b[i] < a[i] && return false
    end
    return false
end

@inline function _lex_less(a::T, b::T) where {T <: xFgh_solution}
    a.sum_violations < b.sum_violations && return true
    b.sum_violations < a.sum_violations && return false
    return _lex_less(a.f, b.f)
end

@inline function _is_dominated_by_front(pop, front, i)
    p = pop[i]
    # most recently added members are closer in the visiting order
    @inbounds for t in length(front):-1:1
        compare(pop[front[t]], p) == 1 && return true
    end
    return false
end

function _non_dominated_sort_ens(pop)
    n = length(pop)
    rank = zeros(Int, n)
    fronts = Vector{Int}[]
    for i in sortperm(pop, lt = _lex_less)
        # being dominated by front k implies being dominated by fronts 1:k-1
        lo, hi = 1, length(fronts) + 1
        while lo < hi
            mid = (lo + hi) >>> 1
            if _is_dominated_by_front(pop, fronts[mid], i)
                lo = mid + 1
            else
                hi = mid
            end
        end
        lo > length(fronts) ? push!(fronts, [i]) : push!(fronts[lo], i)
        rank[i] = lo
    end
    return rank
end

function _non_dominated_sort_naive(pop)
    n = length(pop)

    dom_list = [Int[] for i in 1:n]
    rank = zeros(Int, n)
    dom_count = zeros(Int, n)

    for i in 1:n
        for j in i + 1:n
            comparison = compare(pop[i], pop[j])
            if comparison == 1
                push!(dom_list[i], j)
                dom_count[j] += 1
            elseif comparison == 2
                push!(dom_list[j], i)
                dom_count[i] += 1
            end
        end
    end

    current = findall(iszero, dom_count)
    rank[current] .= 1
    k = 1
    next = Int[]
    while !isempty(current)
        k += 1
        empty!(next)
        for p in current, q in dom_list[p]
            dom_count[q] -= 1
            if dom_count[q] == 0
                rank[q] = k
                push!(next, q)
            end
        end
        current, next = next, current
    end

    return rank
end

"""
    get_fronts(population, computed_ranks = true)

Return each sub-front in an array. If `computed_ranks == true`, this method assumes that
`fast_non_dominated_sort!(population)` has been called before.
"""
function get_fronts(population, computed_ranks = true)
    if computed_ranks
        ranks = [s.rank for s in population]
    else
        ranks = non_dominated_sort(population)
    end

    isempty(ranks) && return Vector{Int}[]
    fronts = [Int[] for i in 1:maximum(ranks)]

    for (i, r) in enumerate(ranks)
        push!(fronts[r], i)
    end

    return fronts
end

"""
    fast_non_dominated_sort!(population)

Sort population using the fast non dominated sorting algorithm. Note that `s.rank` is
updated for each solution `s ∈ population`.
"""
function fast_non_dominated_sort!(pop)
    rank = non_dominated_sort(pop)

    for i in eachindex(pop)
        pop[i].rank = rank[i]
    end

    sort!(pop, by = s -> s.rank, alg = Base.Sort.QuickSort)
    pop
end

"""
    get_non_dominated_solutions_perm(population)

Return a vector of integers `v` such that `population[v]` are the non dominated
solutions contained in `population`.
"""
function get_non_dominated_solutions_perm(population)
    n = length(population)
    n == 0 && return Int[]
    ids = Int[1]

    for i in 2:n
        j = 1
        while j <= length(ids)
            jj = ids[j]
            relation = compare(population[i], population[jj])

            if relation == 2
                # j dominates i
                break
            elseif relation == 1
                # i dominates j
                deleteat!(ids, j)
                continue
            end

            j += 1
        end

        if j > length(ids)
            push!(ids, i)
        end

    end

    return ids
end

"""
    get_non_dominated_solutions(population)

Return the non dominated solutions contained in `population`.
"""
function get_non_dominated_solutions(population)
    return population[get_non_dominated_solutions_perm(population)]
end
